#!/bin/bash
# run_test.sh — SoloVar CDH1 smoke test
#
# Tests the MAF post-processing pipeline on a small synthetic CDH1 MAF:
#   slim_maf → filter_basic → filter_biological
#
# NOTE: filter_af.py is NOT tested here — it requires pysam and
#       bgzipped+indexed VCF files produced by the full annotation pipeline.
#
# Requirements: python3, pandas, numpy  (no pysam needed)
#
# Usage:
#   cd SoloVar/test/
#   bash run_test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLOVAR_DIR="$(dirname "$SCRIPT_DIR")"
SCRIPTS="$SOLOVAR_DIR/scripts"
DATA_DIR="$SCRIPT_DIR/data"
OUT_DIR="$SCRIPT_DIR/output"

mkdir -p "$OUT_DIR"

echo "=== SoloVar CDH1 Smoke Test ==="
echo "Scripts : $SCRIPTS"
echo "Data    : $DATA_DIR"
echo "Output  : $OUT_DIR"
echo ""

# ---- Step 1: slim_maf — coalesce columns, compute VAF ----
echo "[1/3] slim_maf.py — column coalescing & VAF"
python3 "$SCRIPTS/slim_maf.py" \
    --input           "$DATA_DIR/cdh1_test.maf" \
    --output          "$OUT_DIR/slim.tsv" \
    --tumor-types-out "$OUT_DIR/tumor_types_long.tsv" \
    --summary-dir     "$OUT_DIR/summaries/"

# ---- Step 2: filter_basic — FILTER / ML conflict annotation ----
echo "[2/3] filter_basic.py — FILTER/ML conflict annotation"
python3 "$SCRIPTS/filter_basic.py" \
    --in     "$OUT_DIR/slim.tsv" \
    --out    "$OUT_DIR/1F_slim.tsv" \
    --report "$OUT_DIR/1F.report"

# ---- Step 3: filter_biological — biological filtering ----
echo "[3/3] filter_biological.py — LOF / OncoKB / hotspot / pathogenicity filter"
python3 "$SCRIPTS/filter_biological.py" \
    --infile             "$OUT_DIR/1F_slim.tsv" \
    --outfile            "$OUT_DIR/2F_slim.tsv" \
    --outfile_breast_top "$OUT_DIR/2F_breast_top.tsv" \
    --outfile_germline   "$OUT_DIR/germline.tsv"

# ---- Summary ----
echo ""
echo "=== Results ==="
printf "%-30s %s\n" "Input variants:"    "$(tail -n +2 "$DATA_DIR/cdh1_test.maf" | wc -l)"
printf "%-30s %s\n" "After slim_maf:"    "$(tail -n +2 "$OUT_DIR/slim.tsv" | wc -l)"
printf "%-30s %s\n" "After filter_basic:"  "$(tail -n +2 "$OUT_DIR/1F_slim.tsv" | wc -l)"
printf "%-30s %s\n" "After filter_biological:" "$(tail -n +2 "$OUT_DIR/2F_slim.tsv" | wc -l)"
printf "%-30s %s\n" "Breast-dominant:"   "$(tail -n +2 "$OUT_DIR/2F_breast_top.tsv" | wc -l)"
echo ""
echo "Output files written to: $OUT_DIR"
echo "Test complete."
