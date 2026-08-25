#!/usr/bin/env bash
set -e

INPUT_EXCEL="sh2_mutation.xlsx"
SHEET_NAME="Sheet1"
OUTPUT_CSV="sh2_mutation_ESM1b_predictions.csv"
# Change this according to your model set
#MODELS=("esm2_t33_650M_UR50D")
#MODELS=("esm1v_t33_650M_UR90S_1" "esm1v_t33_650M_UR90S_2" "esm1v_t33_650M_UR90S_3" "esm1v_t33_650M_UR90S_4" "esm1v_t33_650M_UR90S_5")
MODELS=("esm1b_t33_650M_UR50S")

python bash.py \
  --input-excel "$INPUT_EXCEL" \
  --sheet-name "$SHEET_NAME" \
  --output-csv "$OUTPUT_CSV" \
  --model-location "${MODELS[@]}" \
  --offset-idx 1