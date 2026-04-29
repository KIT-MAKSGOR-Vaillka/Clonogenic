#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import cv2
import numpy as np


@dataclass(frozen=True)
class Circle:
    x: int
    y: int
    radius: int


@dataclass(frozen=True)
class PhotoConfig:
    min_area: int = 32
    max_area: int = 900
    count_radius_scale: float = 1.0
    hough_param2: float = 38.0
    min_radius_fraction: float = 0.15
    max_radius_fraction: float = 0.29
    min_center_distance_fraction: float = 0.31
    blackhat_quantile: float = 0.982
    blackhat_floor: float = 0.020
    grow_threshold_scale: float = 0.55
    grow_threshold_floor: float = 0.012
    min_seed_darkness: float = 0.055
    min_grow_darkness: float = 0.035
    min_mean_blackhat: float = 0.24
    min_integrated_blackhat: float = 0.0
    strong_dark_gray: float = 0.16
    saturation_floor: float = 0.34
    saturation_delta: float = 0.04
    purple_delta: float = 0.018
    absolute_dark_gray: float = 0.20
    contour_width: int = 2
    clump_stain_saturation_floor: float = 0.42
    clump_stain_saturation_delta: float = 0.07
    clump_stain_purple_delta: float = 0.030
    clump_stain_dark_delta: float = 0.045
    large_clump_max_area_factor: float = 4.0
    split_peak_rel_threshold: float = 0.62
    split_peak_abs_threshold: float = 5.0
    split_peak_min_area: int = 10
    split_max_lobes: int = 4
    split_max_area_per_lobe_factor: float = 1.7


@dataclass
class SegmentRecord:
    image_name: str
    well_index: int
    side: str
    replicate: int
    sample_name: str | None
    concentration: float | None
    dose_gy: float | None
    colony_count: int
    component_kind: str
    colony_label: str
    area_px: int
    centroid_x: float
    centroid_y: float
    bbox_x: int
    bbox_y: int
    bbox_width: int
    bbox_height: int
    circularity: float
    elongation: float
    aspect_ratio: float
    extent: float
    edge_min_px: float
    edge_centroid_px: float
    mean_gray: float
    mean_purple: float
    mean_saturation: float
    mean_blackhat: float
    integrated_blackhat: float


@dataclass
class WellSummary:
    image_name: str
    well_index: int
    side: str
    replicate: int
    sample_name: str | None
    concentration: float | None
    dose_gy: float | None
    counted_colonies: int
    center_x_px: int
    center_y_px: int
    hough_radius_px: int
    count_radius_px: float
    blackhat_threshold: float
    grow_threshold: float
    median_gray: float
    median_purple: float
    median_saturation: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze photographed six-well clonogenic assay plates."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input image paths or glob patterns, for example '*.jpg'.",
    )
    parser.add_argument(
        "--output-dir",
        default="photo_analysis_output",
        help="Directory where output artifacts will be written.",
    )
    parser.add_argument(
        "--min-area",
        type=int,
        default=PhotoConfig.min_area,
        help=(
            "Minimum connected component area in pixels. This is the main size "
            "threshold for moving from all visible specks to countable colonies."
        ),
    )
    parser.add_argument(
        "--max-area",
        type=int,
        default=PhotoConfig.max_area,
        help=(
            "Normal single-colony area limit in pixels. Larger strongly stained "
            "components are counted as merged clumps instead of being discarded."
        ),
    )
    parser.add_argument(
        "--count-radius-scale",
        type=float,
        default=PhotoConfig.count_radius_scale,
        help=(
            "Fraction of the detected well radius used for counting. The default "
            "uses the full detected radius so edge colonies are still eligible."
        ),
    )
    parser.add_argument(
        "--min-mean-blackhat",
        type=float,
        default=PhotoConfig.min_mean_blackhat,
        help=(
            "Minimum mean local-darkness signal for a component. Raise this to "
            "remove weak background/rim marks; lower it if pale colonies disappear."
        ),
    )
    parser.add_argument(
        "--min-integrated-blackhat",
        type=float,
        default=PhotoConfig.min_integrated_blackhat,
        help=(
            "Minimum area multiplied by mean local-darkness signal. This is useful "
            "when both size and intensity should contribute to the cutoff."
        ),
    )
    parser.add_argument(
        "--hough-param2",
        type=float,
        default=PhotoConfig.hough_param2,
        help="OpenCV Hough circle accumulator threshold; lower finds more candidate wells.",
    )
    parser.add_argument(
        "--replicate-layout",
        choices=["auto", "split-columns", "all-wells"],
        default="auto",
        help=(
            "How filename concentrations are assigned to wells. auto uses all six wells "
            "as one condition when the filename has one concentration, and left/right "
            "three-well groups when it has two concentrations."
        ),
    )
    return parser.parse_args()


def expand_inputs(patterns: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        path = Path(pattern)
        if path.exists():
            files.append(path)
        else:
            files.extend(sorted(Path.cwd().glob(pattern)))
    unique = sorted({path.resolve() for path in files})
    if not unique:
        raise SystemExit("No input images matched the provided paths or glob patterns.")
    return [Path(path) for path in unique]


def parse_decimal(raw: str) -> float | None:
    try:
        return float(raw.strip().replace(",", "."))
    except ValueError:
        return None


def split_capture_suffix(stem: str) -> tuple[str, int]:
    match = re.search(r"_(\d+)$", stem)
    if not match:
        return stem, 1
    return stem[: match.start()], int(match.group(1))


def parse_photo_metadata(image_path: Path) -> dict[str, object]:
    plate_key, capture_index = split_capture_suffix(image_path.stem)
    metadata: dict[str, object] = {
        "image_name": image_path.name,
        "plate_key": plate_key,
        "capture_index": capture_index,
        "sample_name": None,
        "dose_gy": None,
        "dose_unit": None,
        "concentrations": [],
        "right_concentration": None,
        "left_concentration": None,
    }

    match = re.search(r"(?i)(^|[_\-\s])([0-9]+(?:[.,][0-9]+)?)\s*gy(?=$|[_\-\s])", plate_key)
    if not match:
        return metadata

    sample_text = plate_key[: match.start()].strip("_- ")
    metadata["sample_name"] = sample_text or None
    metadata["dose_gy"] = parse_decimal(match.group(2))
    metadata["dose_unit"] = "Gy"
    concentration_text = plate_key[match.end() :].strip("_- ")
    concentration_tokens = [token for token in re.split(r"[_\s-]+", concentration_text) if token]
    concentrations = [
        value
        for value in (parse_decimal(token) for token in concentration_tokens)
        if value is not None
    ]
    metadata["concentrations"] = concentrations
    if len(concentration_tokens) >= 1:
        metadata["right_concentration"] = parse_decimal(concentration_tokens[0])
    if len(concentration_tokens) >= 2:
        metadata["left_concentration"] = parse_decimal(concentration_tokens[1])
    return metadata


def well_layout_metadata(
    well_index: int,
    metadata: dict[str, object],
    replicate_layout: str,
) -> tuple[str, int, float | None, float | None]:
    concentrations = metadata.get("concentrations")
    all_wells_layout = replicate_layout == "all-wells" or (
        replicate_layout == "auto" and isinstance(concentrations, list) and len(concentrations) == 1
    )
    if all_wells_layout:
        concentration = concentrations[0] if isinstance(concentrations, list) and concentrations else metadata.get("right_concentration")
        return "all", well_index, concentration, metadata.get("dose_gy")
    if well_index in (1, 3, 5):
        return "left", (well_index + 1) // 2, metadata.get("left_concentration"), metadata.get("dose_gy")
    return "right", well_index // 2, metadata.get("right_concentration"), metadata.get("dose_gy")


def odd_kernel_size(value: float, minimum: int, maximum: int) -> int:
    size = int(round(value))
    size = max(minimum, min(maximum, size))
    if size % 2 == 0:
        size += 1
    return size


def order_circles(circles: list[Circle]) -> list[Circle]:
    if len(circles) != 6:
        raise RuntimeError(f"Expected 6 wells, detected {len(circles)}.")

    x_values = np.asarray([circle.x for circle in circles], dtype=float)
    y_values = np.asarray([circle.y for circle in circles], dtype=float)
    vertical_layout = (y_values.max() - y_values.min()) > (x_values.max() - x_values.min()) * 1.05

    sorted_by_y = sorted(circles, key=lambda circle: circle.y)
    if vertical_layout:
        rows = [sorted(sorted_by_y[index : index + 2], key=lambda circle: circle.x) for index in range(0, 6, 2)]
    else:
        rows = [sorted(sorted_by_y[index : index + 3], key=lambda circle: circle.x) for index in range(0, 6, 3)]
    return [circle for row in rows for circle in row]


def non_overlapping_circles(raw_circles: np.ndarray, shape: tuple[int, int], min_distance: float) -> list[Circle]:
    height, width = shape
    selected: list[Circle] = []
    for raw_x, raw_y, raw_radius in np.round(raw_circles).astype(int):
        circle = Circle(int(raw_x), int(raw_y), int(raw_radius))
        if (
            circle.x - circle.radius * 0.65 < 0
            or circle.y - circle.radius * 0.65 < 0
            or circle.x + circle.radius * 0.65 >= width
            or circle.y + circle.radius * 0.65 >= height
        ):
            continue
        if any(math.hypot(circle.x - current.x, circle.y - current.y) < min_distance for current in selected):
            continue
        selected.append(circle)
        if len(selected) == 6:
            break
    return selected


def detect_wells(image_bgr: np.ndarray, config: PhotoConfig) -> list[Circle]:
    height, width = image_bgr.shape[:2]
    min_dim = min(height, width)
    min_radius = int(round(min_dim * config.min_radius_fraction))
    max_radius = int(round(min_dim * config.max_radius_fraction))
    min_distance = min_dim * config.min_center_distance_fraction

    gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(gray)
    blurred = cv2.GaussianBlur(clahe, (9, 9), 2)

    param2_values = [
        config.hough_param2 + 7,
        config.hough_param2,
        config.hough_param2 - 5,
        config.hough_param2 - 10,
        config.hough_param2 - 15,
    ]

    best: list[Circle] = []
    for param2 in dict.fromkeys(value for value in param2_values if value > 0):
        raw = cv2.HoughCircles(
            blurred,
            cv2.HOUGH_GRADIENT,
            dp=1.2,
            minDist=min_distance,
            param1=80,
            param2=float(param2),
            minRadius=min_radius,
            maxRadius=max_radius,
        )
        if raw is None:
            continue
        selected = non_overlapping_circles(raw[0], (height, width), min_distance * 0.75)
        if len(selected) == 6:
            return order_circles(selected)
        if len(selected) > len(best):
            best = selected

    if len(best) != 6:
        raise RuntimeError(f"Expected 6 wells, detected {len(best)}.")
    return order_circles(best)


def compute_feature_maps(image_bgr: np.ndarray, median_well_radius: float) -> dict[str, np.ndarray]:
    gray_u8 = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)
    gray = gray_u8.astype(np.float32) / 255.0
    rgb = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB).astype(np.float32) / 255.0
    hsv = cv2.cvtColor((rgb * 255).astype(np.uint8), cv2.COLOR_RGB2HSV).astype(np.float32)

    red = rgb[..., 0]
    green = rgb[..., 1]
    blue = rgb[..., 2]
    purple = (red + blue) * 0.5 - green
    saturation = hsv[..., 1] / 255.0

    small_kernel = odd_kernel_size(median_well_radius * 0.075, minimum=9, maximum=51)
    large_kernel = odd_kernel_size(median_well_radius * 0.135, minimum=15, maximum=81)
    small_bh = cv2.morphologyEx(
        gray_u8,
        cv2.MORPH_BLACKHAT,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (small_kernel, small_kernel)),
    ).astype(np.float32) / 255.0
    large_bh = cv2.morphologyEx(
        gray_u8,
        cv2.MORPH_BLACKHAT,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (large_kernel, large_kernel)),
    ).astype(np.float32) / 255.0
    blackhat = cv2.GaussianBlur(np.maximum(small_bh, large_bh), (3, 3), 0.5)

    return {
        "gray": gray,
        "purple": purple,
        "saturation": saturation,
        "blackhat": blackhat,
    }


def component_elongation(xs: np.ndarray, ys: np.ndarray) -> float:
    if xs.size < 3:
        return 1.0
    coords = np.column_stack([xs, ys]).astype(float)
    covariance = np.cov(coords, rowvar=False)
    eigenvalues = np.linalg.eigvalsh(covariance)
    return float(math.sqrt((eigenvalues[-1] + 1e-6) / (eigenvalues[0] + 1e-6)))


def contour_circularity(component: np.ndarray, area_px: int) -> float:
    contours, _ = cv2.findContours(component.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    perimeter = sum(cv2.arcLength(contour, True) for contour in contours)
    if perimeter <= 0:
        return 0.0
    return float(4.0 * math.pi * area_px / (perimeter * perimeter))


def is_component_rejected(
    *,
    area_px: int,
    aspect_ratio: float,
    circularity: float,
    elongation: float,
    extent: float,
    edge_min_px: float,
    circle_radius_px: int,
    mean_gray: float,
    mean_purple: float,
    mean_saturation: float,
    well_median_purple: float,
    mean_blackhat: float,
    integrated_blackhat: float,
    config: PhotoConfig,
) -> bool:
    is_large_component = area_px > config.max_area
    strong_stained = (
        mean_saturation > config.clump_stain_saturation_floor
        and mean_purple > well_median_purple + config.clump_stain_purple_delta * 0.75
        and mean_gray < 0.42
    )
    if mean_blackhat < config.min_mean_blackhat and not (is_large_component and strong_stained):
        return True
    if integrated_blackhat < config.min_integrated_blackhat and not (is_large_component and strong_stained):
        return True
    if is_large_component:
        if not strong_stained:
            return True
        if edge_min_px < circle_radius_px * 0.020 and mean_gray > config.strong_dark_gray:
            if circularity < 0.16 or elongation > 4.8 or aspect_ratio > 4.8:
                return True
        if extent < 0.10 and mean_gray > 0.15:
            return True
        return False
    if aspect_ratio > 5.5 and area_px < 80:
        return True
    if area_px > 80 and elongation > 5.0 and mean_gray > config.strong_dark_gray:
        return True
    if area_px > 130 and (circularity < 0.13 or extent < 0.22 or elongation > 4.0):
        if not (mean_gray < config.strong_dark_gray and mean_saturation > 0.55):
            return True
    if edge_min_px < circle_radius_px * 0.035 and (circularity < 0.24 or elongation > 3.0 or aspect_ratio > 3.0):
        strong_boundary_colony = mean_gray < 0.14 and mean_saturation > 0.55 and area_px < 160
        if not strong_boundary_colony:
            return True
    if edge_min_px < circle_radius_px * 0.015 and mean_purple < well_median_purple + 0.022 and mean_gray > 0.20:
        return True
    return False


def estimate_component_colony_count(
    *,
    component: np.ndarray,
    area_px: int,
    bbox_x: int,
    bbox_y: int,
    bbox_width: int,
    bbox_height: int,
    config: PhotoConfig,
) -> tuple[int, str]:
    if area_px <= config.max_area:
        return 1, "single"

    mega_clump_area = int(round(config.max_area * config.large_clump_max_area_factor))
    if area_px >= mega_clump_area:
        return 1, "large_clump"

    crop = component[bbox_y : bbox_y + bbox_height, bbox_x : bbox_x + bbox_width].astype(np.uint8)
    distance = cv2.distanceTransform(crop, cv2.DIST_L2, 5)
    max_distance = float(distance.max())
    if max_distance < config.split_peak_abs_threshold:
        return 1, "large_clump"

    peak_threshold = max(config.split_peak_abs_threshold, max_distance * config.split_peak_rel_threshold)
    peak_mask = distance >= peak_threshold
    peak_count, _, peak_stats, _ = cv2.connectedComponentsWithStats(peak_mask.astype(np.uint8), 8)
    peaks = [
        label
        for label in range(1, peak_count)
        if int(peak_stats[label, cv2.CC_STAT_AREA]) >= config.split_peak_min_area
    ]
    lobe_count = len(peaks)
    if not (2 <= lobe_count <= config.split_max_lobes):
        return 1, "large_clump"

    area_per_lobe = float(area_px) / float(lobe_count)
    if area_per_lobe < config.min_area * 1.25:
        return 1, "large_clump"
    if area_per_lobe > config.max_area * config.split_max_area_per_lobe_factor:
        return 1, "large_clump"
    return lobe_count, "split_clump"


def assign_colony_labels(records: list[SegmentRecord], well_index: int) -> None:
    next_number = 1
    ordered = sorted(records, key=lambda record: (record.centroid_y, record.centroid_x))
    for record in ordered:
        if record.colony_count <= 1:
            suffix = str(next_number)
        else:
            suffix = f"{next_number}-{next_number + record.colony_count - 1}"
        record.colony_label = f"W{well_index}.{suffix}"
        next_number += max(1, record.colony_count)


def analyze_well(
    *,
    image_name: str,
    well_index: int,
    side: str,
    replicate: int,
    sample_name: str | None,
    concentration: float | None,
    dose_gy: float | None,
    circle: Circle,
    feature_maps: dict[str, np.ndarray],
    xx: np.ndarray,
    yy: np.ndarray,
    config: PhotoConfig,
) -> tuple[WellSummary, list[SegmentRecord], np.ndarray]:
    gray = feature_maps["gray"]
    purple = feature_maps["purple"]
    saturation = feature_maps["saturation"]
    blackhat = feature_maps["blackhat"]

    count_radius = float(circle.radius * config.count_radius_scale)
    distance = np.sqrt((xx - circle.x) ** 2 + (yy - circle.y) ** 2)
    well_mask = distance <= count_radius
    core_mask = distance <= circle.radius * 0.70

    median_gray = float(np.median(gray[core_mask]))
    median_purple = float(np.median(purple[core_mask]))
    median_saturation = float(np.median(saturation[core_mask]))

    dark_abs = np.clip(median_gray - gray, 0.0, None)
    stain_like = (
        (saturation > max(config.saturation_floor, median_saturation + config.saturation_delta))
        | (purple > median_purple + config.purple_delta)
        | (gray < config.absolute_dark_gray)
    )
    clump_stain_mask = (
        (saturation > max(config.clump_stain_saturation_floor, median_saturation + config.clump_stain_saturation_delta))
        & (purple > median_purple + config.clump_stain_purple_delta)
        & (gray < median_gray - config.clump_stain_dark_delta)
        & well_mask
    )
    clump_stain_mask = cv2.morphologyEx(
        clump_stain_mask.astype(np.uint8),
        cv2.MORPH_OPEN,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)),
    ).astype(bool)

    blackhat_threshold = max(float(np.quantile(blackhat[well_mask], config.blackhat_quantile)), config.blackhat_floor)
    grow_threshold = max(blackhat_threshold * config.grow_threshold_scale, config.grow_threshold_floor)

    seed_mask = (
        (blackhat > blackhat_threshold)
        & (dark_abs > config.min_seed_darkness)
        & stain_like
        & well_mask
    ) | clump_stain_mask
    grow_mask = (
        (blackhat > grow_threshold)
        & (dark_abs > config.min_grow_darkness)
        & (stain_like | (gray < median_gray - 0.11))
        & well_mask
    ) | clump_stain_mask

    component_count, labels, stats, centroids = cv2.connectedComponentsWithStats(grow_mask.astype(np.uint8), 8)
    counted_mask = np.zeros(grow_mask.shape, dtype=np.uint8)
    records: list[SegmentRecord] = []

    for label in range(1, component_count):
        area_px = int(stats[label, cv2.CC_STAT_AREA])
        if area_px < config.min_area:
            continue

        component = labels == label
        if not np.any(seed_mask[component]):
            continue

        bbox_x = int(stats[label, cv2.CC_STAT_LEFT])
        bbox_y = int(stats[label, cv2.CC_STAT_TOP])
        bbox_width = int(stats[label, cv2.CC_STAT_WIDTH])
        bbox_height = int(stats[label, cv2.CC_STAT_HEIGHT])
        aspect_ratio = float(max(bbox_width, bbox_height) / max(1, min(bbox_width, bbox_height)))
        extent = float(area_px / max(1, bbox_width * bbox_height))
        circularity = contour_circularity(component, area_px)

        ys, xs = np.where(component)
        elongation = component_elongation(xs, ys)
        distances = np.sqrt((xs - circle.x) ** 2 + (ys - circle.y) ** 2)
        edge_min = float(count_radius - distances.max())
        centroid_x = float(centroids[label][0])
        centroid_y = float(centroids[label][1])
        edge_centroid = float(count_radius - math.hypot(centroid_x - circle.x, centroid_y - circle.y))

        mean_gray = float(gray[component].mean())
        mean_purple = float(purple[component].mean())
        mean_saturation = float(saturation[component].mean())
        mean_blackhat = float(blackhat[component].mean())
        integrated_blackhat = float(area_px * mean_blackhat)

        if is_component_rejected(
            area_px=area_px,
            aspect_ratio=aspect_ratio,
            circularity=circularity,
            elongation=elongation,
            extent=extent,
            edge_min_px=edge_min,
            circle_radius_px=circle.radius,
            mean_gray=mean_gray,
            mean_purple=mean_purple,
            mean_saturation=mean_saturation,
            well_median_purple=median_purple,
            mean_blackhat=mean_blackhat,
            integrated_blackhat=integrated_blackhat,
            config=config,
        ):
            continue

        colony_count, component_kind = estimate_component_colony_count(
            component=component,
            area_px=area_px,
            bbox_x=bbox_x,
            bbox_y=bbox_y,
            bbox_width=bbox_width,
            bbox_height=bbox_height,
            config=config,
        )
        counted_mask[component] = 1
        records.append(
            SegmentRecord(
                image_name=image_name,
                well_index=well_index,
                side=side,
                replicate=replicate,
                sample_name=sample_name,
                concentration=concentration,
                dose_gy=dose_gy,
                colony_count=colony_count,
                component_kind=component_kind,
                colony_label="",
                area_px=area_px,
                centroid_x=round(centroid_x, 3),
                centroid_y=round(centroid_y, 3),
                bbox_x=bbox_x,
                bbox_y=bbox_y,
                bbox_width=bbox_width,
                bbox_height=bbox_height,
                circularity=round(circularity, 6),
                elongation=round(elongation, 6),
                aspect_ratio=round(aspect_ratio, 6),
                extent=round(extent, 6),
                edge_min_px=round(edge_min, 3),
                edge_centroid_px=round(edge_centroid, 3),
                mean_gray=round(mean_gray, 6),
                mean_purple=round(mean_purple, 6),
                mean_saturation=round(mean_saturation, 6),
                mean_blackhat=round(mean_blackhat, 6),
                integrated_blackhat=round(integrated_blackhat, 6),
            )
        )

    assign_colony_labels(records, well_index)
    summary = WellSummary(
        image_name=image_name,
        well_index=well_index,
        side=side,
        replicate=replicate,
        sample_name=sample_name,
        concentration=concentration,
        dose_gy=dose_gy,
        counted_colonies=int(sum(record.colony_count for record in records)),
        center_x_px=circle.x,
        center_y_px=circle.y,
        hough_radius_px=circle.radius,
        count_radius_px=round(count_radius, 3),
        blackhat_threshold=round(blackhat_threshold, 6),
        grow_threshold=round(grow_threshold, 6),
        median_gray=round(median_gray, 6),
        median_purple=round(median_purple, 6),
        median_saturation=round(median_saturation, 6),
    )
    return summary, records, counted_mask


def crop_bounds(circle: Circle, image_shape: tuple[int, int], scale: float = 1.10) -> tuple[int, int, int, int]:
    height, width = image_shape
    radius = int(round(circle.radius * scale))
    left = max(0, circle.x - radius)
    top = max(0, circle.y - radius)
    right = min(width, circle.x + radius)
    bottom = min(height, circle.y + radius)
    return left, top, right, bottom


def save_masked_images(
    *,
    image_bgr: np.ndarray,
    overlay_bgr: np.ndarray,
    all_mask: np.ndarray,
    circles: list[Circle],
    image_output_dir: Path,
) -> None:
    image_output_dir.mkdir(parents=True, exist_ok=True)

    cv2.imwrite(str(image_output_dir / "plate_overlay.png"), overlay_bgr)
    cv2.imwrite(str(image_output_dir / "colony_mask.png"), all_mask * 255)

    colonies_only = np.zeros_like(image_bgr)
    colonies_only[all_mask.astype(bool)] = image_bgr[all_mask.astype(bool)]
    cv2.imwrite(str(image_output_dir / "colonies_only.png"), colonies_only)

    transparent = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2BGRA)
    transparent[..., 3] = all_mask * 255
    cv2.imwrite(str(image_output_dir / "colonies_only_transparent.png"), transparent)

    for index, circle in enumerate(circles, start=1):
        left, top, right, bottom = crop_bounds(circle, image_bgr.shape[:2])
        local_mask = all_mask[top:bottom, left:right]

        cv2.imwrite(str(image_output_dir / f"well_{index:02d}_overlay.png"), overlay_bgr[top:bottom, left:right])
        cv2.imwrite(str(image_output_dir / f"well_{index:02d}_mask.png"), local_mask * 255)

        local_colonies = np.zeros_like(image_bgr[top:bottom, left:right])
        local_colonies[local_mask.astype(bool)] = image_bgr[top:bottom, left:right][local_mask.astype(bool)]
        cv2.imwrite(str(image_output_dir / f"well_{index:02d}_colonies_only.png"), local_colonies)


def draw_overlay(
    image_bgr: np.ndarray,
    circles: list[Circle],
    summaries: list[WellSummary],
    segments: list[SegmentRecord],
    all_mask: np.ndarray,
    config: PhotoConfig,
) -> np.ndarray:
    overlay = image_bgr.copy()
    for circle, summary in zip(circles, summaries, strict=True):
        count_radius = int(round(circle.radius * config.count_radius_scale))
        cv2.circle(overlay, (circle.x, circle.y), count_radius, (255, 220, 0), 2)
        label_x = max(8, int(circle.x - circle.radius * 0.78))
        label_y = max(26, int(circle.y - circle.radius * 0.75))
        cv2.putText(
            overlay,
            f"W{summary.well_index}: {summary.counted_colonies}",
            (label_x, label_y),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.8,
            (0, 255, 255),
            2,
            cv2.LINE_AA,
        )

    contours, _ = cv2.findContours(all_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(overlay, contours, -1, (0, 0, 255), config.contour_width)
    for segment in segments:
        label = segment.colony_label.split(".", 1)[-1]
        anchor = (int(round(segment.centroid_x)) + 4, int(round(segment.centroid_y)) - 4)
        cv2.putText(
            overlay,
            label,
            anchor,
            cv2.FONT_HERSHEY_SIMPLEX,
            0.38,
            (0, 0, 0),
            2,
            cv2.LINE_AA,
        )
        cv2.putText(
            overlay,
            label,
            anchor,
            cv2.FONT_HERSHEY_SIMPLEX,
            0.38,
            (255, 255, 0),
            1,
            cv2.LINE_AA,
        )
    return overlay


def write_csv(path: Path, rows: list[object]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(asdict(rows[0]).keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_dict_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_group_summary(well_summaries: list[WellSummary], output_path: Path) -> None:
    groups: dict[float | None, list[WellSummary]] = {}
    for summary in well_summaries:
        groups.setdefault(summary.concentration, []).append(summary)

    rows: list[dict[str, object]] = []
    for concentration, summaries in groups.items():
        counts = [int(summary.counted_colonies) for summary in summaries]
        values = np.asarray(counts, dtype=float)
        row_names = ",".join(sorted({summary.side for summary in summaries}))
        rows.append(
            {
                "row_name": row_names,
                "concentration": concentration,
                "replicates": len(counts),
                "mean_count": round(float(values.mean()), 3),
                "std_count": round(float(values.std(ddof=0)), 3),
                "counts": ",".join(str(value) for value in counts),
            }
        )
    rows.sort(key=lambda row: str(row["concentration"]))
    write_dict_csv(output_path, rows)


def analyze_image(
    image_path: Path,
    output_root: Path,
    config: PhotoConfig,
    replicate_layout: str,
) -> tuple[list[WellSummary], list[SegmentRecord]]:
    image_bgr = cv2.imread(str(image_path))
    if image_bgr is None:
        raise RuntimeError(f"Could not read image: {image_path}")

    metadata = parse_photo_metadata(image_path)
    circles = detect_wells(image_bgr, config)
    median_radius = float(np.median([circle.radius for circle in circles]))
    feature_maps = compute_feature_maps(image_bgr, median_radius)

    height, width = image_bgr.shape[:2]
    yy, xx = np.indices((height, width))
    image_output_dir = output_root / image_path.stem
    image_output_dir.mkdir(parents=True, exist_ok=True)

    summaries: list[WellSummary] = []
    all_segments: list[SegmentRecord] = []
    all_mask = np.zeros((height, width), dtype=np.uint8)

    for well_index, circle in enumerate(circles, start=1):
        side, replicate, concentration, dose_gy = well_layout_metadata(well_index, metadata, replicate_layout)
        summary, segments, local_mask = analyze_well(
            image_name=image_path.name,
            well_index=well_index,
            side=side,
            replicate=replicate,
            sample_name=metadata.get("sample_name") if isinstance(metadata.get("sample_name"), str) else None,
            concentration=concentration if isinstance(concentration, float) else None,
            dose_gy=dose_gy if isinstance(dose_gy, float) else None,
            circle=circle,
            feature_maps=feature_maps,
            xx=xx,
            yy=yy,
            config=config,
        )
        summaries.append(summary)
        all_segments.extend(segments)
        all_mask[local_mask.astype(bool)] = 1

    overlay = draw_overlay(image_bgr, circles, summaries, all_segments, all_mask, config)
    save_masked_images(
        image_bgr=image_bgr,
        overlay_bgr=overlay,
        all_mask=all_mask,
        circles=circles,
        image_output_dir=image_output_dir,
    )

    write_csv(image_output_dir / "well_counts.csv", summaries)
    write_csv(image_output_dir / "segments.csv", all_segments)
    write_group_summary(summaries, image_output_dir / "group_summary.csv")
    with (image_output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "image": str(image_path),
                "metadata": metadata,
                "config": asdict(config),
                "total_counted_colonies": int(sum(summary.counted_colonies for summary in summaries)),
                "wells": [asdict(summary) for summary in summaries],
            },
            handle,
            ensure_ascii=False,
            indent=2,
        )

    return summaries, all_segments


def main() -> None:
    args = parse_args()
    config = PhotoConfig(
        min_area=args.min_area,
        max_area=args.max_area,
        count_radius_scale=args.count_radius_scale,
        hough_param2=args.hough_param2,
        min_mean_blackhat=args.min_mean_blackhat,
        min_integrated_blackhat=args.min_integrated_blackhat,
    )
    output_root = Path(args.output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    all_summaries: list[WellSummary] = []
    all_segments: list[SegmentRecord] = []
    for image_path in expand_inputs(args.inputs):
        summaries, segments = analyze_image(image_path, output_root, config, args.replicate_layout)
        all_summaries.extend(summaries)
        all_segments.extend(segments)

    write_csv(output_root / "batch_well_counts.csv", all_summaries)
    write_csv(output_root / "batch_segments.csv", all_segments)


if __name__ == "__main__":
    main()
