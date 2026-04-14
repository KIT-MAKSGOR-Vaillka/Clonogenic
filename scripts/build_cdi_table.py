#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
ET.register_namespace("", NS)
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TEMPLATE = SCRIPT_DIR.parent / "templates" / "CDI UiO-66 (Zr).xlsx"


@dataclass
class MatrixBlock:
    cell_line: str
    doses: list[float]
    concentrations: list[float]
    raw: dict[float, dict[float, float]]
    normalized: dict[float, dict[float, float]]
    cdi_raw: dict[float, dict[float, float]]
    cdi_100: dict[float, dict[float, float]]
    material: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build CDI verification tables and machine-readable CSV from clonogenic TIFF scan outputs."
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="LABEL=DIR",
        help="Cell-line label and analysis root, for example 4T1=analysis_output_java_scans_final",
    )
    parser.add_argument(
        "--material",
        default="Material",
        help="Material label to write into the long-form CSV and use in heatmap titles.",
    )
    parser.add_argument(
        "--template",
        default=str(DEFAULT_TEMPLATE),
        help="XLSX template to preserve formatting for the verification workbook.",
    )
    parser.add_argument(
        "--output-prefix",
        default="cdi_report",
        help="Prefix for generated files.",
    )
    return parser.parse_args()


def parse_input_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"Expected LABEL=DIR, got: {spec}")
    label, raw_dir = spec.split("=", 1)
    label = label.strip()
    path = Path(raw_dir.strip()).resolve()
    if not label:
        raise ValueError(f"Missing label in input spec: {spec}")
    if not path.exists():
        raise FileNotFoundError(path)
    return label, path


def sort_concentrations(values: set[float]) -> list[float]:
    zero = [value for value in values if abs(value) < 1e-9]
    nonzero = sorted((value for value in values if abs(value) >= 1e-9), reverse=True)
    return zero + nonzero


def sort_nonzero(values: set[float]) -> list[float]:
    return sorted((value for value in values if abs(value) >= 1e-9), reverse=True)


def format_number(value: float | None, decimals: int = 3) -> str:
    if value is None:
        return ""
    text = f"{value:.{decimals}f}"
    text = text.rstrip("0").rstrip(".")
    return text or "0"


def format_header_dose(value: float) -> str:
    return f"{format_number(value, 3)} Gy"


def load_scan_root(label: str, root: Path, material: str) -> MatrixBlock:
    summary_paths = sorted(path for path in root.glob("*/summary.json") if path.is_file())
    if not summary_paths:
        raise FileNotFoundError(f"No summary.json files found under {root}")

    values: dict[float, dict[float, float]] = {}
    doses: set[float] = set()
    concentrations: set[float] = set()

    for summary_path in summary_paths:
        with summary_path.open(encoding="utf-8") as handle:
            summary = json.load(handle)
        dose = float(summary["metadata"]["dose_gy"])
        doses.add(dose)

        group_path = summary_path.parent / "group_summary.csv"
        with group_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                concentration = float(row["concentration"])
                mean_count = float(row["mean_count"])
                concentrations.add(concentration)
                values.setdefault(concentration, {})[dose] = mean_count

    dose_order = sorted(doses)
    concentration_order = sort_concentrations(concentrations)
    control = values[0.0][0.0]

    normalized: dict[float, dict[float, float]] = {}
    cdi_raw: dict[float, dict[float, float]] = {}
    cdi_100: dict[float, dict[float, float]] = {}

    for concentration in concentration_order:
        normalized[concentration] = {}
        for dose in dose_order:
            raw_value = values.get(concentration, {}).get(dose)
            if raw_value is None:
                continue
            normalized[concentration][dose] = raw_value / control * 100.0

    for concentration in concentration_order:
        if abs(concentration) < 1e-9:
            continue
        cdi_raw[concentration] = {}
        cdi_100[concentration] = {}
        drug_alone_raw = values[concentration][0.0]
        drug_alone_norm = normalized[concentration][0.0]
        for dose in dose_order:
            if abs(dose) < 1e-9:
                continue
            radiation_alone_raw = values[0.0][dose]
            combined_raw = values[concentration][dose]
            radiation_alone_norm = normalized[0.0][dose]
            combined_norm = normalized[concentration][dose]
            cdi_raw[concentration][dose] = combined_raw / (drug_alone_raw * radiation_alone_raw)
            cdi_100[concentration][dose] = combined_norm * 100.0 / (drug_alone_norm * radiation_alone_norm)

    return MatrixBlock(
        cell_line=label,
        doses=dose_order,
        concentrations=concentration_order,
        raw=values,
        normalized=normalized,
        cdi_raw=cdi_raw,
        cdi_100=cdi_100,
        material=material,
    )


def write_long_csv(blocks: list[MatrixBlock], output_path: Path) -> None:
    fieldnames = [
        "material",
        "cell_line",
        "dose_gy",
        "concentration_mg_ml",
        "raw_mean_count",
        "normalized_pct",
        "cdi_raw",
        "cdi_100",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for block in blocks:
            for concentration in block.concentrations:
                for dose in block.doses:
                    writer.writerow(
                        {
                            "material": block.material,
                            "cell_line": block.cell_line,
                            "dose_gy": format_number(dose, 3),
                            "concentration_mg_ml": format_number(concentration, 3),
                            "raw_mean_count": format_number(block.raw.get(concentration, {}).get(dose), 6),
                            "normalized_pct": format_number(block.normalized.get(concentration, {}).get(dose), 6),
                            "cdi_raw": format_number(block.cdi_raw.get(concentration, {}).get(dose), 12),
                            "cdi_100": format_number(block.cdi_100.get(concentration, {}).get(dose), 12),
                        }
                    )


def build_table_assignments(blocks: list[MatrixBlock]) -> OrderedDict[str, str | float | None]:
    if len(blocks) > 2:
        raise ValueError("The verification table layout supports at most two cell lines per sheet.")

    block_specs = [
        {
            "raw_label_col": "A",
            "raw_data_cols": ["B", "C", "D", "E"],
            "norm_label_col": "P",
            "norm_data_cols": ["Q", "R", "S", "T"],
        },
        {
            "raw_label_col": "H",
            "raw_data_cols": ["I", "J", "K", "L"],
            "norm_label_col": "W",
            "norm_data_cols": ["X", "Y", "Z", "AA"],
        },
    ]

    assignments: OrderedDict[str, str | float | None] = OrderedDict()

    for spec in block_specs:
        label_col = spec["raw_label_col"]
        norm_label_col = spec["norm_label_col"]
        for ref in [f"{label_col}1", f"{norm_label_col}1", f"{label_col}9", f"{label_col}13", f"{norm_label_col}9", f"{norm_label_col}13"]:
            assignments[ref] = ""
        for row in range(3, 7):
            assignments[f"{label_col}{row}"] = ""
            assignments[f"{norm_label_col}{row}"] = ""
        for col in spec["raw_data_cols"]:
            assignments[f"{col}2"] = ""
            for row in range(3, 7):
                assignments[f"{col}{row}"] = ""
            assignments[f"{col}8"] = ""
            for row in range(9, 12):
                assignments[f"{col}{row}"] = ""
            for row in range(13, 16):
                assignments[f"{col}{row}"] = ""
        for col in spec["norm_data_cols"]:
            assignments[f"{col}2"] = ""
            for row in range(3, 7):
                assignments[f"{col}{row}"] = ""
            assignments[f"{col}8"] = ""
            for row in range(9, 12):
                assignments[f"{col}{row}"] = ""
            for row in range(13, 16):
                assignments[f"{col}{row}"] = ""

    for index, block in enumerate(blocks):
        spec = block_specs[index]
        raw_label_col = spec["raw_label_col"]
        raw_data_cols = spec["raw_data_cols"]
        norm_label_col = spec["norm_label_col"]
        norm_data_cols = spec["norm_data_cols"]

        raw_concs = block.concentrations
        cdi_concs = sort_nonzero(set(block.concentrations))
        nonzero_doses = [dose for dose in block.doses if abs(dose) >= 1e-9]

        assignments[f"{raw_label_col}1"] = block.cell_line
        assignments[f"{norm_label_col}1"] = block.cell_line
        assignments[f"{raw_label_col}9"] = "CDI"
        assignments[f"{raw_label_col}13"] = "CDI (100)"
        assignments[f"{norm_label_col}9"] = "CDI"
        assignments[f"{norm_label_col}13"] = "CDI (100)"

        for row_index, concentration in enumerate(raw_concs, start=3):
            assignments[f"{raw_label_col}{row_index}"] = concentration
            assignments[f"{norm_label_col}{row_index}"] = concentration

        for col, dose in zip(raw_data_cols, block.doses):
            assignments[f"{col}2"] = format_header_dose(dose)
        for col, dose in zip(norm_data_cols, block.doses):
            assignments[f"{col}2"] = format_header_dose(dose)

        for row_index, concentration in enumerate(raw_concs, start=3):
            for col, dose in zip(raw_data_cols, block.doses):
                assignments[f"{col}{row_index}"] = block.raw.get(concentration, {}).get(dose)
            for col, dose in zip(norm_data_cols, block.doses):
                assignments[f"{col}{row_index}"] = block.normalized.get(concentration, {}).get(dose)

        for col, dose in zip(raw_data_cols, nonzero_doses):
            assignments[f"{col}8"] = format_header_dose(dose)
        for col, dose in zip(norm_data_cols, nonzero_doses):
            assignments[f"{col}8"] = format_header_dose(dose)

        for row_offset, concentration in enumerate(cdi_concs):
            row_raw = 9 + row_offset
            row_norm = 13 + row_offset
            for col, dose in zip(raw_data_cols, nonzero_doses):
                assignments[f"{col}{row_raw}"] = block.cdi_raw.get(concentration, {}).get(dose)
                assignments[f"{col}{row_norm}"] = block.cdi_100.get(concentration, {}).get(dose)
            for col, dose in zip(norm_data_cols, nonzero_doses):
                assignments[f"{col}{row_raw}"] = block.cdi_raw.get(concentration, {}).get(dose)
                assignments[f"{col}{row_norm}"] = block.cdi_100.get(concentration, {}).get(dose)

    return assignments


def write_verification_csv(assignments: OrderedDict[str, str | float | None], output_path: Path) -> None:
    refs = list(assignments.keys())
    max_row = max(int(re.search(r"(\d+)$", ref).group(1)) for ref in refs)
    max_col = max(column_to_index(re.search(r"([A-Z]+)", ref).group(1)) for ref in refs)

    grid = [["" for _ in range(max_col)] for _ in range(max_row)]
    for ref, value in assignments.items():
        row_index = int(re.search(r"(\d+)$", ref).group(1)) - 1
        col_index = column_to_index(re.search(r"([A-Z]+)", ref).group(1)) - 1
        if value in ("", None):
            grid[row_index][col_index] = ""
        elif isinstance(value, str):
            grid[row_index][col_index] = value
        else:
            grid[row_index][col_index] = format_number(float(value), 12)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(grid)


def column_to_index(label: str) -> int:
    value = 0
    for char in label:
        value = value * 26 + (ord(char) - ord("A") + 1)
    return value


def find_or_create_cell(sheet_root: ET.Element, ref: str) -> ET.Element:
    row_number = int(re.search(r"(\d+)$", ref).group(1))
    cell_column = re.search(r"([A-Z]+)", ref).group(1)
    sheet_data = sheet_root.find(f"{{{NS}}}sheetData")
    row_element = None
    for candidate in sheet_data.findall(f"{{{NS}}}row"):
        if int(candidate.attrib["r"]) == row_number:
            row_element = candidate
            break
    if row_element is None:
        row_element = ET.SubElement(sheet_data, f"{{{NS}}}row", {"r": str(row_number)})

    for candidate in row_element.findall(f"{{{NS}}}c"):
        if candidate.attrib.get("r") == ref:
            return candidate

    cell = ET.SubElement(row_element, f"{{{NS}}}c", {"r": ref})
    cells = row_element.findall(f"{{{NS}}}c")
    cells.sort(key=lambda item: column_to_index(re.search(r"([A-Z]+)", item.attrib["r"]).group(1)))
    for child in list(row_element):
        row_element.remove(child)
    for child in cells:
        row_element.append(child)
    return cell


def set_cell_value(cell: ET.Element, value: str | float | None) -> None:
    for child in list(cell):
        cell.remove(child)
    if value in ("", None):
        cell.attrib.pop("t", None)
        return
    if isinstance(value, str):
        cell.attrib["t"] = "inlineStr"
        is_node = ET.SubElement(cell, f"{{{NS}}}is")
        text_node = ET.SubElement(is_node, f"{{{NS}}}t")
        text_node.text = value
    else:
        cell.attrib.pop("t", None)
        value_node = ET.SubElement(cell, f"{{{NS}}}v")
        value_node.text = format_number(float(value), 12)


def write_verification_xlsx(template_path: Path, assignments: OrderedDict[str, str | float | None], output_path: Path) -> None:
    with zipfile.ZipFile(template_path, "r") as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
        sheet_root = ET.fromstring(members["xl/worksheets/sheet1.xml"])
        for ref, value in assignments.items():
            cell = find_or_create_cell(sheet_root, ref)
            set_cell_value(cell, value)
        members["xl/worksheets/sheet1.xml"] = ET.tostring(sheet_root, encoding="utf-8", xml_declaration=True)

    with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        with zipfile.ZipFile(temp_path, "w") as archive:
            for name, data in members.items():
                archive.writestr(name, data)
        shutil.move(temp_path, output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    inputs = [parse_input_spec(spec) for spec in args.input]
    blocks = [load_scan_root(label, path, args.material) for label, path in inputs]

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    long_csv_path = output_prefix.with_name(output_prefix.name + "_long.csv")
    verification_csv_path = output_prefix.with_name(output_prefix.name + "_verification.csv")
    verification_xlsx_path = output_prefix.with_name(output_prefix.name + "_verification.xlsx")

    write_long_csv(blocks, long_csv_path)
    assignments = build_table_assignments(blocks)
    write_verification_csv(assignments, verification_csv_path)

    template_path = Path(args.template)
    if template_path.exists():
        write_verification_xlsx(template_path, assignments, verification_xlsx_path)
    else:
        print(f"Template not found: {template_path}. Skipped XLSX output.")

    print(f"Wrote {long_csv_path}")
    print(f"Wrote {verification_csv_path}")
    if template_path.exists():
        print(f"Wrote {verification_xlsx_path}")


if __name__ == "__main__":
    main()
