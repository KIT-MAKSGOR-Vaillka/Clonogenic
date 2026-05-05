#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import asdict, dataclass, replace
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
    profile: str = "universal"
    min_area: int = 32
    max_area: int = 900
    count_radius_scale: float = 1.0
    edge_rescue_radius_scale: float = 1.06
    reference_well_radius_px: float = 240.0
    small_well_scale_start_fraction: float = 0.90
    small_well_area_factor: float = 0.75
    min_scaled_area: int = 8
    low_resolution_min_area: int = 11
    low_resolution_min_core_radius: float = 3.4
    small_noise_area_factor: float = 2.1
    small_noise_integrated_blackhat: float = 16.0
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
    min_mean_blackhat: float = 0.20
    min_integrated_blackhat: float = 0.0
    strong_dark_gray: float = 0.16
    saturation_floor: float = 0.34
    saturation_delta: float = 0.04
    purple_delta: float = 0.018
    absolute_dark_gray: float = 0.20
    contour_width: int = 1
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
    core_radius_px: float
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
    min_area_px: int
    max_area_px: int
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
        "--edge-rescue-radius-scale",
        type=float,
        default=PhotoConfig.edge_rescue_radius_scale,
        help=(
            "Expanded search radius used only for low-saturation photo wells. "
            "This rescues stained colonies on the physical rim without changing "
            "normal high-saturation images."
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
        "--small-noise-integrated-blackhat",
        type=float,
        default=PhotoConfig.small_noise_integrated_blackhat,
        help=(
            "Small-component integrated darkness cutoff used in noisy/low-saturation "
            "photo wells. Raise to remove tiny specks; lower if real small colonies disappear."
        ),
    )
    parser.add_argument(
        "--contour-width",
        type=int,
        default=PhotoConfig.contour_width,
        help="Contour line width for colony overlays.",
    )
    parser.add_argument(
        "--hough-param2",
        type=float,
        default=PhotoConfig.hough_param2,
        help="OpenCV Hough circle accumulator threshold; lower finds more candidate wells.",
    )
    parser.add_argument(
        "--profile",
        choices=["universal", "auto", "default", "zr"],
        default=PhotoConfig.profile,
        help=(
            "Compatibility option. Detection is now universal and does not branch "
            "by sample name; legacy values auto/default/zr are accepted."
        ),
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
    capture_index = int(match.group(1))
    if capture_index <= 1:
        return stem, 1
    return stem[: match.start()], capture_index


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


def config_for_metadata(config: PhotoConfig, metadata: dict[str, object]) -> PhotoConfig:
    del metadata
    return replace(config, profile="universal")


def scaled_area_threshold(
    area_px: int,
    circle_radius_px: int,
    config: PhotoConfig,
    *,
    floor_px: int,
    scale_up: bool = False,
) -> int:
    scale = (float(circle_radius_px) / config.reference_well_radius_px) ** 2
    if not scale_up:
        if circle_radius_px >= config.reference_well_radius_px * config.small_well_scale_start_fraction:
            scale = 1.0
        else:
            scale = min(1.0, scale)
            scale *= config.small_well_area_factor
    return max(floor_px, int(round(area_px * scale)))


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


def component_core_radius(
    component: np.ndarray,
    bbox_x: int,
    bbox_y: int,
    bbox_width: int,
    bbox_height: int,
) -> float:
    crop = component[bbox_y : bbox_y + bbox_height, bbox_x : bbox_x + bbox_width].astype(np.uint8)
    if crop.size == 0:
        return 0.0
    distance = cv2.distanceTransform(crop, cv2.DIST_L2, 5)
    return float(distance.max())


def low_extent_component_mask(mask: np.ndarray, *, min_area_px: int, max_extent: float) -> np.ndarray:
    component_count, labels, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8), 8)
    rejected = np.zeros(mask.shape, dtype=bool)
    for label in range(1, component_count):
        area_px = int(stats[label, cv2.CC_STAT_AREA])
        if area_px < min_area_px:
            continue
        bbox_width = int(stats[label, cv2.CC_STAT_WIDTH])
        bbox_height = int(stats[label, cv2.CC_STAT_HEIGHT])
        extent = float(area_px / max(1, bbox_width * bbox_height))
        if extent < max_extent:
            rejected[labels == label] = True
    return rejected


def is_component_rejected(
    *,
    area_px: int,
    min_area_px: int,
    max_area_px: int,
    aspect_ratio: float,
    circularity: float,
    elongation: float,
    extent: float,
    core_radius_px: float,
    edge_min_px: float,
    circle_radius_px: int,
    mean_gray: float,
    mean_purple: float,
    mean_saturation: float,
    well_median_purple: float,
    well_median_saturation: float,
    mean_blackhat: float,
    integrated_blackhat: float,
    config: PhotoConfig,
) -> bool:
    is_large_component = area_px > max_area_px
    color_signal = (
        mean_saturation > max(config.saturation_floor, well_median_saturation + config.saturation_delta)
        or mean_purple > well_median_purple + config.purple_delta
        or mean_gray < config.absolute_dark_gray
    )
    stain_color_signal = (
        mean_saturation > max(config.saturation_floor, well_median_saturation + config.saturation_delta)
        or (
            mean_saturation > max(0.26, well_median_saturation + 0.035)
            and mean_purple > well_median_purple + config.purple_delta
        )
        or mean_purple > well_median_purple + config.purple_delta * 1.65
    )
    solid_dark_colony = (
        area_px >= min_area_px
        and mean_gray < 0.23
        and mean_saturation > max(0.32, well_median_saturation + 0.055)
        and mean_purple > well_median_purple + 0.014
        and circularity > 0.32
        and elongation < 2.9
        and extent > 0.38
        and core_radius_px >= max(2.2, circle_radius_px * 0.007)
    )
    strong_edge_blob = (
        mean_gray < 0.18
        and mean_saturation > max(0.55, well_median_saturation + 0.18)
        and mean_purple > well_median_purple + 0.040
        and core_radius_px >= max(6.0, circle_radius_px * 0.020)
        and extent > 0.25
        and elongation < 3.0
    )
    low_saturation_background = well_median_saturation < 0.24
    low_resolution_well = min_area_px < config.min_area
    needs_noise_rescue = low_saturation_background or low_resolution_well
    needs_small_noise_cleanup = needs_noise_rescue and not low_resolution_well
    low_res_compact_colony = False
    low_res_borderline_compact_colony = False
    low_res_pale_compact_colony = False
    if low_resolution_well:
        low_res_compact_colony = (
            mean_gray < 0.34
            and mean_blackhat > config.min_mean_blackhat * 1.08
            and mean_saturation > well_median_saturation + 0.095
            and (
                mean_saturation > well_median_saturation + 0.13
                or mean_purple > well_median_purple + 0.002
            )
            and circularity > 0.32
            and elongation < 3.6
            and extent > 0.36
        )
        low_res_borderline_compact_colony = (
            low_res_compact_colony
            and mean_gray < 0.32
            and mean_blackhat > config.min_mean_blackhat * 1.18
            and circularity > 0.50
            and extent > 0.52
        )
        low_res_pale_compact_colony = (
            area_px >= max(18, int(round(config.low_resolution_min_area * 1.6)))
            and area_px <= 180
            and mean_gray < 0.40
            and mean_blackhat > config.min_mean_blackhat * 0.85
            and (
                mean_saturation > well_median_saturation + 0.045
                or mean_purple > well_median_purple + 0.0015
            )
            and circularity > 0.22
            and elongation < 3.8
            and extent > 0.28
            and core_radius_px >= 1.3
        )
        if area_px < config.low_resolution_min_area:
            low_res_borderline_area = area_px >= max(10, config.low_resolution_min_area - 1)
            if not (low_res_borderline_area and low_res_borderline_compact_colony):
                return True
        if area_px < int(round(config.low_resolution_min_area * 2.2)):
            if core_radius_px < 2.8 and not low_res_compact_colony:
                return True
    weak_low_chroma = (
        mean_gray > 0.34
        and mean_saturation < max(0.34, well_median_saturation + 0.03)
        and mean_purple < well_median_purple + 0.022
    )
    small_noise_area_limit = max(
        int(round(min_area_px * config.small_noise_area_factor)),
        48 if min_area_px >= config.min_area else 22,
    )
    small_noise_integrated_floor = scaled_area_threshold(
        int(round(config.small_noise_integrated_blackhat)),
        circle_radius_px,
        config,
        floor_px=6,
        scale_up=False,
    )
    if (
        needs_small_noise_cleanup
        and area_px < small_noise_area_limit
        and integrated_blackhat < small_noise_integrated_floor
        and not solid_dark_colony
    ):
        return True
    if needs_small_noise_cleanup and area_px < max(72, int(round(min_area_px * 2.25))):
        strong_tiny_colony = (
            mean_gray < 0.15
            and mean_saturation > max(0.58, well_median_saturation + 0.22)
            and mean_purple > well_median_purple + 0.025
            and integrated_blackhat >= 20.0
            and circularity > 0.48
            and extent > 0.52
            and core_radius_px >= 3.4
        )
        if not strong_tiny_colony:
            return True
    if needs_small_noise_cleanup and weak_low_chroma and area_px < max(80, int(round(max_area_px * 0.35))):
        return True
    if needs_small_noise_cleanup and area_px < max(int(round(min_area_px * 1.8)), 18):
        compact_colony = (
            color_signal
            and circularity > 0.28
            and elongation < 2.6
            and (mean_blackhat > config.min_mean_blackhat * 0.65 or solid_dark_colony)
        )
        if not compact_colony:
            return True
    if needs_small_noise_cleanup and area_px < max(int(round(min_area_px * 3.0)), 34) and mean_gray > 0.30 and mean_saturation < 0.34:
        return True

    edge_compact_colony = (
        low_saturation_background
        and edge_min_px < circle_radius_px * 0.10
        and area_px <= max(240, int(round(max_area_px * 0.35)))
        and mean_gray < 0.40
        and mean_blackhat > config.min_mean_blackhat * 0.70
        and (
            mean_saturation > max(0.18, well_median_saturation + 0.055)
            or mean_purple > max(0.012, well_median_purple + 0.006)
            or mean_gray < 0.26
        )
        and circularity > 0.16
        and elongation < 4.3
        and aspect_ratio < 4.0
        and extent > 0.24
        and core_radius_px >= 1.3
    )

    if low_saturation_background:
        thin_rim_arc = (
            edge_min_px < circle_radius_px * 0.16
            and area_px > max(96, int(round(min_area_px * 3.0)))
            and core_radius_px < max(3.4, circle_radius_px * 0.012)
            and (circularity < 0.34 or elongation > 2.5 or extent < 0.28)
        )
        if thin_rim_arc and not edge_compact_colony:
            return True
        rim_like_large = (
            edge_min_px < circle_radius_px * 0.13
            and area_px > max(140, int(round(max_area_px * 0.55)))
            and (circularity < 0.20 or extent < 0.22 or elongation > 4.0 or aspect_ratio > 3.5)
        )
    else:
        rim_like_large = (
            edge_min_px < circle_radius_px * 0.055
            and area_px > 220
            and (circularity < 0.26 or extent < 0.18 or elongation > 3.2 or aspect_ratio > 3.4)
        )
    if rim_like_large and not strong_edge_blob and not edge_compact_colony:
        return True

    low_chroma_large_artifact = (
        low_saturation_background
        and area_px > max(120, int(round(min_area_px * 4.0)))
        and mean_gray > 0.24
        and mean_saturation < max(0.20, well_median_saturation + 0.060)
        and mean_purple < max(0.010, well_median_purple + 0.008)
        and not strong_edge_blob
        and (circularity < 0.40 or extent < 0.46 or elongation > 1.9 or core_radius_px < 4.5)
    )
    if low_chroma_large_artifact:
        return True

    strong_stained = (
        mean_saturation > config.clump_stain_saturation_floor
        and mean_purple > well_median_purple + config.clump_stain_purple_delta * 0.75
        and mean_gray < 0.42
    )
    if (
        mean_blackhat < config.min_mean_blackhat
        and not (is_large_component and strong_stained)
        and not low_res_pale_compact_colony
        and not edge_compact_colony
    ):
        return True
    if integrated_blackhat < config.min_integrated_blackhat and not (is_large_component and strong_stained):
        return True
    if (
        needs_noise_rescue
        and not stain_color_signal
        and not low_res_pale_compact_colony
        and not edge_compact_colony
        and area_px < max(270, int(round(max_area_px * 0.30)))
    ):
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
        if low_saturation_background:
            strong_boundary_colony = (
                edge_compact_colony
                or (
                    color_signal
                    and area_px < max(180, int(round(max_area_px * 0.70)))
                    and extent > 0.20
                    and elongation < 3.6
                    and (circularity > 0.18 or mean_saturation > 0.48 or mean_gray < 0.16)
                )
            )
        else:
            strong_boundary_colony = mean_gray < 0.14 and mean_saturation > 0.55 and area_px < 160
        if not strong_boundary_colony:
            return True
    if (
        edge_min_px < circle_radius_px * 0.015
        and mean_purple < well_median_purple + 0.022
        and mean_gray > 0.20
        and not edge_compact_colony
    ):
        return True
    return False


def is_edge_subcomponent_rescue(
    *,
    area_px: int,
    min_area_px: int,
    max_area_px: int,
    aspect_ratio: float,
    circularity: float,
    elongation: float,
    extent: float,
    core_radius_px: float,
    edge_min_px: float,
    circle_radius_px: int,
    mean_gray: float,
    mean_purple: float,
    mean_saturation: float,
    well_median_purple: float,
    well_median_saturation: float,
    mean_blackhat: float,
    config: PhotoConfig,
) -> bool:
    if area_px < max(min_area_px, 18):
        return False
    if area_px > max(int(round(max_area_px * 1.45)), 1300):
        return False
    if edge_min_px > circle_radius_px * 0.10:
        return False
    small_edge_subcomponent = area_px < 50
    if small_edge_subcomponent:
        if core_radius_px < 1.3:
            return False
        if extent < 0.26 or elongation > 4.2 or aspect_ratio > 3.8:
            return False
        if circularity < 0.16:
            return False
    elif core_radius_px < max(2.4, circle_radius_px * 0.008):
        return False
    if extent < 0.14 or elongation > 5.2 or aspect_ratio > 4.6:
        return False
    if circularity < 0.075 and core_radius_px < max(4.2, circle_radius_px * 0.014):
        return False
    rescue_color = (
        mean_saturation > well_median_saturation + 0.060
        or mean_purple > well_median_purple + 0.004
    )
    if mean_gray > 0.40 and not rescue_color:
        return False
    if circularity < 0.13 and (elongation > 3.2 or extent < 0.28) and mean_saturation < 0.45:
        return False
    edge_color_or_signal = (
        rescue_color
        or mean_blackhat > config.min_mean_blackhat * 0.55
    )
    return mean_gray < 0.44 and edge_color_or_signal


def estimate_component_colony_count(
    *,
    component: np.ndarray,
    area_px: int,
    min_area_px: int,
    max_area_px: int,
    bbox_x: int,
    bbox_y: int,
    bbox_width: int,
    bbox_height: int,
    config: PhotoConfig,
) -> tuple[int, str]:
    if area_px <= max_area_px:
        return 1, "single"

    mega_clump_area = int(round(max_area_px * config.large_clump_max_area_factor))
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
    if area_per_lobe < min_area_px * 1.25:
        return 1, "large_clump"
    if area_per_lobe > max_area_px * config.split_max_area_per_lobe_factor:
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

    min_area_px = scaled_area_threshold(
        config.min_area,
        circle.radius,
        config,
        floor_px=config.min_scaled_area,
        scale_up=False,
    )
    high_dose_size_cleanup = False
    max_area_px = config.max_area
    distance = np.sqrt((xx - circle.x) ** 2 + (yy - circle.y) ** 2)
    base_count_radius = float(circle.radius * config.count_radius_scale)
    core_mask = distance <= circle.radius * 0.70

    median_gray = float(np.median(gray[core_mask]))
    median_purple = float(np.median(purple[core_mask]))
    median_saturation = float(np.median(saturation[core_mask]))
    low_saturation_strength = float(np.clip((0.34 - median_saturation) / 0.12, 0.0, 1.0))
    low_saturation_background = low_saturation_strength >= 0.35
    detection_config = replace(
        config,
        min_mean_blackhat=0.24 - ((0.24 - config.min_mean_blackhat) * low_saturation_strength),
    )
    high_dose_size_cleanup = (
        dose_gy is not None
        and dose_gy >= 4.0
        and low_saturation_background
    )
    if high_dose_size_cleanup:
        min_area_px = max(
            min_area_px,
            scaled_area_threshold(
                36,
                circle.radius,
                config,
                floor_px=16,
                scale_up=False,
            ),
        )
    count_radius = base_count_radius
    if low_saturation_background:
        count_radius = max(count_radius, float(circle.radius * config.edge_rescue_radius_scale))
    threshold_mask = distance <= base_count_radius
    well_mask = distance <= count_radius

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
    solid_colony_seed = np.zeros(gray.shape, dtype=bool)
    solid_colony_grow = np.zeros(gray.shape, dtype=bool)
    if low_saturation_background:
        solid_colony_seed = (
            (gray < min(0.23, median_gray - 0.18))
            & (saturation > max(0.34, median_saturation + 0.070))
            & (purple > median_purple + 0.014)
            & well_mask
        )
        solid_colony_grow = (
            (gray < min(0.29, median_gray - 0.12))
            & (saturation > max(0.28, median_saturation + 0.040))
            & (purple > median_purple + 0.010)
            & well_mask
        )

    blackhat_threshold = max(float(np.quantile(blackhat[threshold_mask], config.blackhat_quantile)), config.blackhat_floor)
    grow_threshold = max(blackhat_threshold * config.grow_threshold_scale, config.grow_threshold_floor)
    edge_band = (
        (distance >= base_count_radius * 0.82)
        & (distance <= base_count_radius + circle.radius * 0.025)
        & well_mask
    )
    edge_dark_seed = (
        low_saturation_background
        & edge_band
        & (blackhat > grow_threshold * 0.70)
        & (dark_abs > 0.14)
        & (gray < min(0.45, median_gray - 0.09))
        & (
            (saturation > median_saturation + 0.030)
            | (purple > median_purple + 0.006)
            | (blackhat > blackhat_threshold)
        )
    )
    rim_seed_artifact = low_extent_component_mask(
        edge_dark_seed,
        min_area_px=max(450, int(round(circle.radius * 1.35))),
        max_extent=0.08,
    )
    if np.any(rim_seed_artifact):
        edge_dark_seed = edge_dark_seed & ~rim_seed_artifact
        rim_break_mask = cv2.dilate(
            rim_seed_artifact.astype(np.uint8),
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)),
        ).astype(bool)
    else:
        rim_break_mask = rim_seed_artifact
    seed_mask = (
        (blackhat > blackhat_threshold)
        & (dark_abs > config.min_seed_darkness)
        & stain_like
        & well_mask
    ) | clump_stain_mask | solid_colony_seed | edge_dark_seed
    grow_mask = (
        (blackhat > grow_threshold)
        & (dark_abs > config.min_grow_darkness)
        & (stain_like | (gray < median_gray - 0.11))
        & well_mask
    ) | clump_stain_mask | solid_colony_grow | edge_dark_seed
    if np.any(rim_break_mask):
        grow_mask = grow_mask & ~rim_break_mask

    component_count, labels, stats, centroids = cv2.connectedComponentsWithStats(grow_mask.astype(np.uint8), 8)
    counted_mask = np.zeros(grow_mask.shape, dtype=np.uint8)
    records: list[SegmentRecord] = []

    for label in range(1, component_count):
        area_px = int(stats[label, cv2.CC_STAT_AREA])
        if area_px < min_area_px:
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
        core_radius = component_core_radius(component, bbox_x, bbox_y, bbox_width, bbox_height)

        ys, xs = np.where(component)
        elongation = component_elongation(xs, ys)
        distances = np.sqrt((xs - circle.x) ** 2 + (ys - circle.y) ** 2)
        edge_min = float(base_count_radius - distances.max())
        centroid_x = float(centroids[label][0])
        centroid_y = float(centroids[label][1])
        edge_centroid = float(base_count_radius - math.hypot(centroid_x - circle.x, centroid_y - circle.y))

        mean_gray = float(gray[component].mean())
        mean_purple = float(purple[component].mean())
        mean_saturation = float(saturation[component].mean())
        mean_blackhat = float(blackhat[component].mean())
        integrated_blackhat = float(area_px * mean_blackhat)

        rejected = is_component_rejected(
            area_px=area_px,
            min_area_px=min_area_px,
            max_area_px=max_area_px,
            aspect_ratio=aspect_ratio,
            circularity=circularity,
            elongation=elongation,
            extent=extent,
            core_radius_px=core_radius,
            edge_min_px=edge_min,
            circle_radius_px=circle.radius,
            mean_gray=mean_gray,
            mean_purple=mean_purple,
            mean_saturation=mean_saturation,
            well_median_purple=median_purple,
            well_median_saturation=median_saturation,
            mean_blackhat=mean_blackhat,
            integrated_blackhat=integrated_blackhat,
            config=detection_config,
        )
        if rejected:
            edge_split_candidate = (
                median_saturation < 0.24
                and area_px > max(500, int(round(max_area_px * 0.70)))
                and edge_min < circle.radius * 0.085
                and extent < 0.18
            )
            if edge_split_candidate:
                sub_source = component & (edge_dark_seed | solid_colony_seed)
                sub_count, sub_labels, sub_stats, sub_centroids = cv2.connectedComponentsWithStats(
                    sub_source.astype(np.uint8),
                    8,
                )
                for sub_label in range(1, sub_count):
                    sub_area_px = int(sub_stats[sub_label, cv2.CC_STAT_AREA])
                    if sub_area_px < max(min_area_px, 50):
                        continue
                    subcomponent = sub_labels == sub_label
                    sub_bbox_x = int(sub_stats[sub_label, cv2.CC_STAT_LEFT])
                    sub_bbox_y = int(sub_stats[sub_label, cv2.CC_STAT_TOP])
                    sub_bbox_width = int(sub_stats[sub_label, cv2.CC_STAT_WIDTH])
                    sub_bbox_height = int(sub_stats[sub_label, cv2.CC_STAT_HEIGHT])
                    sub_aspect_ratio = float(
                        max(sub_bbox_width, sub_bbox_height) / max(1, min(sub_bbox_width, sub_bbox_height))
                    )
                    sub_extent = float(sub_area_px / max(1, sub_bbox_width * sub_bbox_height))
                    sub_circularity = contour_circularity(subcomponent, sub_area_px)
                    sub_core_radius = component_core_radius(
                        subcomponent,
                        sub_bbox_x,
                        sub_bbox_y,
                        sub_bbox_width,
                        sub_bbox_height,
                    )
                    sub_ys, sub_xs = np.where(subcomponent)
                    sub_elongation = component_elongation(sub_xs, sub_ys)
                    sub_distances = np.sqrt((sub_xs - circle.x) ** 2 + (sub_ys - circle.y) ** 2)
                    sub_edge_min = float(base_count_radius - sub_distances.max())
                    sub_centroid_x = float(sub_centroids[sub_label][0])
                    sub_centroid_y = float(sub_centroids[sub_label][1])
                    sub_edge_centroid = float(
                        base_count_radius - math.hypot(sub_centroid_x - circle.x, sub_centroid_y - circle.y)
                    )
                    sub_mean_gray = float(gray[subcomponent].mean())
                    sub_mean_purple = float(purple[subcomponent].mean())
                    sub_mean_saturation = float(saturation[subcomponent].mean())
                    sub_mean_blackhat = float(blackhat[subcomponent].mean())
                    sub_integrated_blackhat = float(sub_area_px * sub_mean_blackhat)
                    if not is_edge_subcomponent_rescue(
                        area_px=sub_area_px,
                        min_area_px=min_area_px,
                        max_area_px=max_area_px,
                        aspect_ratio=sub_aspect_ratio,
                        circularity=sub_circularity,
                        elongation=sub_elongation,
                        extent=sub_extent,
                        core_radius_px=sub_core_radius,
                        edge_min_px=sub_edge_min,
                        circle_radius_px=circle.radius,
                        mean_gray=sub_mean_gray,
                        mean_purple=sub_mean_purple,
                        mean_saturation=sub_mean_saturation,
                        well_median_purple=median_purple,
                        well_median_saturation=median_saturation,
                        mean_blackhat=sub_mean_blackhat,
                            config=detection_config,
                    ):
                        continue
                    counted_mask[subcomponent] = 1
                    records.append(
                        SegmentRecord(
                            image_name=image_name,
                            well_index=well_index,
                            side=side,
                            replicate=replicate,
                            sample_name=sample_name,
                            concentration=concentration,
                            dose_gy=dose_gy,
                            colony_count=1,
                            component_kind="edge_rescue",
                            colony_label="",
                            area_px=sub_area_px,
                            centroid_x=round(sub_centroid_x, 3),
                            centroid_y=round(sub_centroid_y, 3),
                            bbox_x=sub_bbox_x,
                            bbox_y=sub_bbox_y,
                            bbox_width=sub_bbox_width,
                            bbox_height=sub_bbox_height,
                            circularity=round(sub_circularity, 6),
                            elongation=round(sub_elongation, 6),
                            aspect_ratio=round(sub_aspect_ratio, 6),
                            extent=round(sub_extent, 6),
                            core_radius_px=round(sub_core_radius, 3),
                            edge_min_px=round(sub_edge_min, 3),
                            edge_centroid_px=round(sub_edge_centroid, 3),
                            mean_gray=round(sub_mean_gray, 6),
                            mean_purple=round(sub_mean_purple, 6),
                            mean_saturation=round(sub_mean_saturation, 6),
                            mean_blackhat=round(sub_mean_blackhat, 6),
                            integrated_blackhat=round(sub_integrated_blackhat, 6),
                        )
                    )
            continue

        colony_count, component_kind = estimate_component_colony_count(
            component=component,
            area_px=area_px,
            min_area_px=min_area_px,
            max_area_px=max_area_px,
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
                core_radius_px=round(core_radius, 3),
                edge_min_px=round(edge_min, 3),
                edge_centroid_px=round(edge_centroid, 3),
                mean_gray=round(mean_gray, 6),
                mean_purple=round(mean_purple, 6),
                mean_saturation=round(mean_saturation, 6),
                mean_blackhat=round(mean_blackhat, 6),
                integrated_blackhat=round(integrated_blackhat, 6),
            )
        )

    if median_saturation < 0.24:
        counted_guard = cv2.dilate(
            counted_mask.astype(np.uint8),
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)),
        ).astype(bool)
        edge_rescue_ring = (
            (distance >= base_count_radius * 0.80)
            & (distance <= count_radius)
            & well_mask
        )
        edge_candidate_mask = (
            edge_rescue_ring
            & ~counted_guard
            & (blackhat > grow_threshold * 0.52)
            & (dark_abs > 0.070)
            & (gray < min(0.44, median_gray - 0.055))
            & (
                (
                    (saturation > max(0.33, median_saturation + 0.100))
                    & (
                        (purple > max(0.014, median_purple + 0.010))
                        | (gray < 0.28)
                    )
                )
                | (
                    (purple > max(0.026, median_purple + 0.018))
                    & (saturation > max(0.22, median_saturation + 0.050))
                )
                | (gray < 0.24)
            )
        )
        edge_candidate_mask = edge_candidate_mask & ~low_extent_component_mask(
            edge_candidate_mask,
            min_area_px=max(80, int(round(circle.radius * 0.50))),
            max_extent=0.12,
        )
        edge_count, edge_labels, edge_stats, edge_centroids = cv2.connectedComponentsWithStats(
            edge_candidate_mask.astype(np.uint8),
            8,
        )
        for edge_label in range(1, edge_count):
            edge_area_px = int(edge_stats[edge_label, cv2.CC_STAT_AREA])
            if edge_area_px < max(min_area_px, 10) or edge_area_px > 320:
                continue
            edge_component = edge_labels == edge_label
            edge_bbox_x = int(edge_stats[edge_label, cv2.CC_STAT_LEFT])
            edge_bbox_y = int(edge_stats[edge_label, cv2.CC_STAT_TOP])
            edge_bbox_width = int(edge_stats[edge_label, cv2.CC_STAT_WIDTH])
            edge_bbox_height = int(edge_stats[edge_label, cv2.CC_STAT_HEIGHT])
            edge_aspect_ratio = float(
                max(edge_bbox_width, edge_bbox_height) / max(1, min(edge_bbox_width, edge_bbox_height))
            )
            edge_extent = float(edge_area_px / max(1, edge_bbox_width * edge_bbox_height))
            edge_circularity = contour_circularity(edge_component, edge_area_px)
            edge_core_radius = component_core_radius(
                edge_component,
                edge_bbox_x,
                edge_bbox_y,
                edge_bbox_width,
                edge_bbox_height,
            )
            edge_ys, edge_xs = np.where(edge_component)
            edge_elongation = component_elongation(edge_xs, edge_ys)
            edge_distances = np.sqrt((edge_xs - circle.x) ** 2 + (edge_ys - circle.y) ** 2)
            edge_min = float(base_count_radius - edge_distances.max())
            if edge_min > circle.radius * 0.10 or edge_min < -circle.radius * 0.08:
                continue
            edge_centroid_x = float(edge_centroids[edge_label][0])
            edge_centroid_y = float(edge_centroids[edge_label][1])
            edge_edge_centroid = float(
                base_count_radius - math.hypot(edge_centroid_x - circle.x, edge_centroid_y - circle.y)
            )
            edge_mean_gray = float(gray[edge_component].mean())
            edge_mean_purple = float(purple[edge_component].mean())
            edge_mean_saturation = float(saturation[edge_component].mean())
            edge_mean_blackhat = float(blackhat[edge_component].mean())
            edge_integrated_blackhat = float(edge_area_px * edge_mean_blackhat)
            if edge_area_px < 60 and (
                edge_circularity < 0.42
                or edge_extent < 0.45
                or edge_elongation > 2.8
                or edge_core_radius < 2.1
            ):
                continue
            edge_color_signal = (
                edge_mean_saturation > max(0.34, median_saturation + 0.100)
                or (
                    edge_mean_purple > max(0.026, median_purple + 0.018)
                    and edge_mean_saturation > max(0.22, median_saturation + 0.050)
                )
                or edge_mean_gray < 0.21
            )
            if edge_mean_gray > 0.34 and not (
                edge_mean_saturation > 0.40
                and edge_mean_purple > max(0.026, median_purple + 0.018)
            ):
                continue
            if edge_min < 0:
                outside_edge_signal = (
                    edge_mean_saturation > 0.42
                    or edge_mean_purple > max(0.030, median_purple + 0.024)
                    or edge_mean_gray < 0.20
                )
                if not outside_edge_signal:
                    continue
                outside_edge_shape = (
                    (
                        edge_circularity > 0.52
                        and edge_extent > 0.50
                        and edge_elongation < 2.8
                    )
                    or (
                        edge_area_px >= 90
                        and edge_core_radius >= 3.4
                        and edge_extent > 0.45
                        and edge_elongation < 2.8
                    )
                )
                if not outside_edge_shape:
                    continue
            edge_rescue_shape = (
                edge_mean_gray < 0.42
                and edge_mean_blackhat > config.min_mean_blackhat * 0.55
                and edge_color_signal
                and (
                    edge_mean_saturation > max(0.30, median_saturation + 0.090)
                    or edge_mean_purple > max(0.016, median_purple + 0.012)
                    or edge_mean_gray < 0.24
                )
                and edge_circularity > 0.12
                and edge_elongation < 4.5
                and edge_aspect_ratio < 4.2
                and edge_extent > 0.20
                and edge_core_radius >= 1.2
            )
            if edge_area_px < 18:
                edge_rescue_shape = (
                    edge_rescue_shape
                    and edge_mean_gray < 0.32
                    and edge_mean_blackhat > config.min_mean_blackhat * 0.90
                    and edge_circularity > 0.45
                    and edge_extent > 0.52
                    and edge_elongation < 2.8
                    and edge_aspect_ratio < 2.5
                )
            if edge_area_px > 120:
                edge_rescue_shape = (
                    edge_rescue_shape
                    and (
                        edge_mean_saturation > 0.40
                        or edge_mean_purple > max(0.034, median_purple + 0.024)
                        or edge_mean_gray < 0.21
                    )
                    and edge_circularity > 0.22
                    and edge_extent > 0.32
                    and edge_elongation < 3.2
                )
            if not edge_rescue_shape:
                continue
            counted_mask[edge_component] = 1
            records.append(
                SegmentRecord(
                    image_name=image_name,
                    well_index=well_index,
                    side=side,
                    replicate=replicate,
                    sample_name=sample_name,
                    concentration=concentration,
                    dose_gy=dose_gy,
                    colony_count=1,
                    component_kind="edge_second_pass",
                    colony_label="",
                    area_px=edge_area_px,
                    centroid_x=round(edge_centroid_x, 3),
                    centroid_y=round(edge_centroid_y, 3),
                    bbox_x=edge_bbox_x,
                    bbox_y=edge_bbox_y,
                    bbox_width=edge_bbox_width,
                    bbox_height=edge_bbox_height,
                    circularity=round(edge_circularity, 6),
                    elongation=round(edge_elongation, 6),
                    aspect_ratio=round(edge_aspect_ratio, 6),
                    extent=round(edge_extent, 6),
                    core_radius_px=round(edge_core_radius, 3),
                    edge_min_px=round(edge_min, 3),
                    edge_centroid_px=round(edge_edge_centroid, 3),
                    mean_gray=round(edge_mean_gray, 6),
                    mean_purple=round(edge_mean_purple, 6),
                    mean_saturation=round(edge_mean_saturation, 6),
                    mean_blackhat=round(edge_mean_blackhat, 6),
                    integrated_blackhat=round(edge_integrated_blackhat, 6),
                )
            )

    if median_saturation < 0.24:
        counted_guard = cv2.dilate(
            counted_mask.astype(np.uint8),
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)),
        ).astype(bool)
        stain_rescue_mask = (
            well_mask
            & ~counted_guard
            & (gray < min(0.44, median_gray - 0.065))
            & (
                (
                    (saturation > max(0.30, median_saturation + 0.085))
                    & (
                        (purple > max(0.012, median_purple + 0.008))
                        | (gray < 0.28)
                    )
                )
                | (
                    (purple > max(0.030, median_purple + 0.022))
                    & (saturation > max(0.20, median_saturation + 0.045))
                )
                | (
                    (gray < 0.16)
                    & (saturation > max(0.26, median_saturation + 0.055))
                )
            )
        )
        stain_rescue_mask = cv2.morphologyEx(
            stain_rescue_mask.astype(np.uint8),
            cv2.MORPH_OPEN,
            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)),
        ).astype(bool)
        stain_rescue_mask = stain_rescue_mask & ~low_extent_component_mask(
            stain_rescue_mask,
            min_area_px=max(120, int(round(circle.radius * 0.40))),
            max_extent=0.12,
        )
        rescue_count, rescue_labels, rescue_stats, rescue_centroids = cv2.connectedComponentsWithStats(
            stain_rescue_mask.astype(np.uint8),
            8,
        )
        min_stain_rescue_area_px = max(32, min_area_px)
        if high_dose_size_cleanup:
            min_stain_rescue_area_px = max(
                min_stain_rescue_area_px,
                scaled_area_threshold(
                    60,
                    circle.radius,
                    config,
                    floor_px=15,
                    scale_up=False,
                ),
            )
        for rescue_label in range(1, rescue_count):
            rescue_area_px = int(rescue_stats[rescue_label, cv2.CC_STAT_AREA])
            if rescue_area_px < min_stain_rescue_area_px or rescue_area_px > 2800:
                continue
            rescue_component = rescue_labels == rescue_label
            rescue_bbox_x = int(rescue_stats[rescue_label, cv2.CC_STAT_LEFT])
            rescue_bbox_y = int(rescue_stats[rescue_label, cv2.CC_STAT_TOP])
            rescue_bbox_width = int(rescue_stats[rescue_label, cv2.CC_STAT_WIDTH])
            rescue_bbox_height = int(rescue_stats[rescue_label, cv2.CC_STAT_HEIGHT])
            rescue_aspect_ratio = float(
                max(rescue_bbox_width, rescue_bbox_height)
                / max(1, min(rescue_bbox_width, rescue_bbox_height))
            )
            rescue_extent = float(rescue_area_px / max(1, rescue_bbox_width * rescue_bbox_height))
            rescue_circularity = contour_circularity(rescue_component, rescue_area_px)
            rescue_core_radius = component_core_radius(
                rescue_component,
                rescue_bbox_x,
                rescue_bbox_y,
                rescue_bbox_width,
                rescue_bbox_height,
            )
            rescue_ys, rescue_xs = np.where(rescue_component)
            rescue_elongation = component_elongation(rescue_xs, rescue_ys)
            rescue_distances = np.sqrt((rescue_xs - circle.x) ** 2 + (rescue_ys - circle.y) ** 2)
            rescue_edge_min = float(base_count_radius - rescue_distances.max())
            if rescue_edge_min < -circle.radius * 0.10:
                continue
            rescue_centroid_x = float(rescue_centroids[rescue_label][0])
            rescue_centroid_y = float(rescue_centroids[rescue_label][1])
            rescue_edge_centroid = float(
                base_count_radius - math.hypot(rescue_centroid_x - circle.x, rescue_centroid_y - circle.y)
            )
            rescue_mean_gray = float(gray[rescue_component].mean())
            rescue_mean_purple = float(purple[rescue_component].mean())
            rescue_mean_saturation = float(saturation[rescue_component].mean())
            rescue_mean_blackhat = float(blackhat[rescue_component].mean())
            rescue_integrated_blackhat = float(rescue_area_px * rescue_mean_blackhat)
            rescue_color = (
                rescue_mean_saturation > max(0.38, median_saturation + 0.15)
                or rescue_mean_purple > max(0.045, median_purple + 0.035)
                or rescue_mean_gray < 0.15
            )
            compact_residual_colony = (
                rescue_color
                and rescue_area_px >= min_stain_rescue_area_px
                and rescue_mean_gray < 0.25
                and rescue_extent > 0.42
                and rescue_elongation < 3.0
                and rescue_aspect_ratio < 3.0
                and rescue_core_radius >= 2.2
                and rescue_circularity > 0.22
            )
            large_residual_colony = (
                rescue_area_px >= 90
                and rescue_mean_gray < 0.28
                and (
                    rescue_mean_saturation > max(0.35, median_saturation + 0.12)
                    or rescue_mean_purple > max(0.040, median_purple + 0.030)
                )
                and rescue_extent > 0.25
                and rescue_elongation < 4.2
                and rescue_aspect_ratio < 4.0
                and rescue_core_radius >= 3.0
                and (rescue_circularity > 0.14 or rescue_area_px >= 220)
            )
            huge_residual_colony = (
                rescue_area_px >= 450
                and rescue_mean_gray < 0.22
                and rescue_mean_saturation > max(0.45, median_saturation + 0.18)
                and rescue_mean_purple > max(0.055, median_purple + 0.045)
                and rescue_extent > 0.22
                and rescue_elongation < 4.6
                and rescue_aspect_ratio < 4.4
                and rescue_core_radius >= 5.0
            )
            if rescue_edge_min < 0 and not (
                huge_residual_colony
                or (
                    rescue_mean_saturation > 0.48
                    and rescue_mean_purple > max(0.050, median_purple + 0.040)
                    and rescue_extent > 0.32
                    and rescue_elongation < 3.4
                )
            ):
                continue
            if not (compact_residual_colony or large_residual_colony or huge_residual_colony):
                continue
            counted_mask[rescue_component] = 1
            records.append(
                SegmentRecord(
                    image_name=image_name,
                    well_index=well_index,
                    side=side,
                    replicate=replicate,
                    sample_name=sample_name,
                    concentration=concentration,
                    dose_gy=dose_gy,
                    colony_count=1,
                    component_kind="stain_rescue",
                    colony_label="",
                    area_px=rescue_area_px,
                    centroid_x=round(rescue_centroid_x, 3),
                    centroid_y=round(rescue_centroid_y, 3),
                    bbox_x=rescue_bbox_x,
                    bbox_y=rescue_bbox_y,
                    bbox_width=rescue_bbox_width,
                    bbox_height=rescue_bbox_height,
                    circularity=round(rescue_circularity, 6),
                    elongation=round(rescue_elongation, 6),
                    aspect_ratio=round(rescue_aspect_ratio, 6),
                    extent=round(rescue_extent, 6),
                    core_radius_px=round(rescue_core_radius, 3),
                    edge_min_px=round(rescue_edge_min, 3),
                    edge_centroid_px=round(rescue_edge_centroid, 3),
                    mean_gray=round(rescue_mean_gray, 6),
                    mean_purple=round(rescue_mean_purple, 6),
                    mean_saturation=round(rescue_mean_saturation, 6),
                    mean_blackhat=round(rescue_mean_blackhat, 6),
                    integrated_blackhat=round(rescue_integrated_blackhat, 6),
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
        min_area_px=min_area_px,
        max_area_px=max_area_px,
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
        count_radius = int(round(summary.count_radius_px))
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
            0.30,
            (0, 0, 0),
            1,
            cv2.LINE_AA,
        )
        cv2.putText(
            overlay,
            label,
            anchor,
            cv2.FONT_HERSHEY_SIMPLEX,
            0.30,
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
    config = config_for_metadata(config, metadata)
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
        profile=args.profile,
        min_area=args.min_area,
        max_area=args.max_area,
        count_radius_scale=args.count_radius_scale,
        edge_rescue_radius_scale=args.edge_rescue_radius_scale,
        hough_param2=args.hough_param2,
        min_mean_blackhat=args.min_mean_blackhat,
        min_integrated_blackhat=args.min_integrated_blackhat,
        small_noise_integrated_blackhat=args.small_noise_integrated_blackhat,
        contour_width=args.contour_width,
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
