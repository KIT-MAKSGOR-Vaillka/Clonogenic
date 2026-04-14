#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  scripts/run_material_pipeline.sh \
    --input-dir DIR \
    --output-dir DIR \
    --material "UiO-66 (Hf)" \
    [--cell-line 4T1] \
    [--template templates/CDI UiO-66 (Zr).xlsx]

Description:
  1. Finds all TIFF scans in --input-dir
  2. Runs the Java scan analyzer
  3. Builds CDI CSV/XLSX verification tables
  4. Renders the heatmap
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
MATERIAL=""
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

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$MATERIAL" ]]; then
  usage >&2
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

mapfile -d '' TIFF_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -iname '*.tif' -o -iname '*.tiff' \) -print0 | sort -z)

if [[ ${#TIFF_FILES[@]} -eq 0 ]]; then
  echo "No TIFF scans found in: $INPUT_DIR" >&2
  exit 1
fi

"$ROOT_DIR/scripts/run_clonogenic_java.sh" --output-dir "$OUTPUT_DIR" "${TIFF_FILES[@]}"

python3 "$ROOT_DIR/scripts/build_cdi_table.py" \
  --input "$CELL_LINE=$OUTPUT_DIR" \
  --material "$MATERIAL" \
  --template "$TEMPLATE_PATH" \
  --output-prefix "$OUTPUT_DIR/cdi_${CELL_LINE,,}"

python3 "$ROOT_DIR/scripts/heatmap.py" \
  --input "$OUTPUT_DIR/cdi_${CELL_LINE,,}_long.csv" \
  --output "$OUTPUT_DIR/cdi_${CELL_LINE,,}_heatmap.png" \
  --material "$MATERIAL"
