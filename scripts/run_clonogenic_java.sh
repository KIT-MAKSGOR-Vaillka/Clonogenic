#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_DIR="$ROOT_DIR/src/java"
BUILD_DIR="$ROOT_DIR/build/classes"

IJ_JAR_DEFAULT="/home/maksegr/Applications/Fiji.app/jars/ij-1.54p.jar"
IJ_JAR="${FIJI_IJ_JAR:-${IJ_JAR:-$IJ_JAR_DEFAULT}}"

if [[ ! -f "$IJ_JAR" ]]; then
  echo "ImageJ/Fiji jar not found: $IJ_JAR" >&2
  echo "Set FIJI_IJ_JAR=/absolute/path/to/ij-*.jar before running this script." >&2
  exit 1
fi

mkdir -p "$BUILD_DIR"

javac \
  -d "$BUILD_DIR" \
  -cp "$IJ_JAR:$SRC_DIR" \
  "$SRC_DIR/WellMaskRefiner.java" \
  "$SRC_DIR/ClonogenicAnalyzer.java"

java -cp "$BUILD_DIR:$IJ_JAR" ClonogenicAnalyzer "$@"
