#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  scripts/run_material_pipeline.sh \
    --input-dir DIR \
    --output-dir DIR \
    [--material "UiO-66 (Hf)"] \
    [--dataset auto|hf|zr] \
    [--config settings.properties] \
    [--layout well_layout.csv] \
    [--manifest plate_manifest.csv] \
    [--skip-cdi] \
    [--cell-line 4T1] \
    [--template templates/CDI UiO-66 (Zr).xlsx]

Description:
  1. Finds all TIFF scans in --input-dir
  2. Runs the Java scan analyzer
  3. Optionally builds CDI CSV/XLSX verification tables
  4. Optionally renders the heatmap
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
MATERIAL=""
DATASET="auto"
CONFIG_PATH=""
LAYOUT_PATH=""
MANIFEST_PATH=""
SKIP_CDI=0
CELL_LINE="4T1"
TEMPLATE_PATH="$ROOT_DIR/templates/CDI UiO-66 (Zr).xlsx"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --material)
      MATERIAL="$2"
      shift 2
      ;;
    --dataset)
      DATASET="$2"
      shift 2
      ;;
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --layout)
      LAYOUT_PATH="$2"
      shift 2
      ;;
    --manifest)
      MANIFEST_PATH="$2"
      shift 2
      ;;
    --skip-cdi)
      SKIP_CDI=1
      shift
      ;;
    --cell-line)
      CELL_LINE="$2"
      shift 2
      ;;
    --template)
      TEMPLATE_PATH="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
  usage >&2
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

case "$DATASET" in
  auto|hf|zr)
    ;;
  *)
    echo "Unsupported --dataset value: $DATASET" >&2
    echo "Expected one of: auto, hf, zr" >&2
    exit 1
    ;;
esac

if [[ -n "$CONFIG_PATH" && ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

if [[ -n "$LAYOUT_PATH" && ! -f "$LAYOUT_PATH" ]]; then
  echo "Layout CSV not found: $LAYOUT_PATH" >&2
  exit 1
fi

if [[ -n "$MANIFEST_PATH" && ! -f "$MANIFEST_PATH" ]]; then
  echo "Manifest CSV not found: $MANIFEST_PATH" >&2
  exit 1
fi

mapfile -d '' TIFF_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -iname '*.tif' -o -iname '*.tiff' \) -print0 | sort -z)

if [[ ${#TIFF_FILES[@]} -eq 0 ]]; then
  echo "No TIFF scans found in: $INPUT_DIR" >&2
  exit 1
fi

JAVA_ARGS=(--output-dir "$OUTPUT_DIR" --dataset "$DATASET")

if [[ -n "$CONFIG_PATH" ]]; then
  JAVA_ARGS+=(--config "$CONFIG_PATH")
fi

if [[ -n "$LAYOUT_PATH" ]]; then
  JAVA_ARGS+=(--layout "$LAYOUT_PATH")
fi

if [[ -n "$MANIFEST_PATH" ]]; then
  JAVA_ARGS+=(--manifest "$MANIFEST_PATH")
fi

JAVA_ARGS+=("${TIFF_FILES[@]}")

"$ROOT_DIR/scripts/run_clonogenic_java.sh" "${JAVA_ARGS[@]}"

if [[ "$SKIP_CDI" -eq 1 ]]; then
  exit 0
fi

if [[ -z "$MATERIAL" ]]; then
  echo "--material is required unless --skip-cdi is used." >&2
  exit 1
fi

python3 "$ROOT_DIR/scripts/build_cdi_table.py" \
  --input "$CELL_LINE=$OUTPUT_DIR" \
  --material "$MATERIAL" \
  --template "$TEMPLATE_PATH" \
  --output-prefix "$OUTPUT_DIR/cdi_${CELL_LINE,,}"

python3 "$ROOT_DIR/scripts/heatmap.py" \
  --input "$OUTPUT_DIR/cdi_${CELL_LINE,,}_long.csv" \
  --output "$OUTPUT_DIR/cdi_${CELL_LINE,,}_heatmap.png" \
  --material "$MATERIAL"
