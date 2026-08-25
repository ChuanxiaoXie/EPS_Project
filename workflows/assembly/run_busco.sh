#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 --input FILE --lineage NAME_OR_DIR --mode genome|proteins|transcriptome --output DIR --run-name NAME [--threads INT] [--offline]" >&2
}

input=""
lineage=""
mode=""
output=""
run_name=""
threads=15
offline=false
while (($#)); do
  case "$1" in
    --input) input=$2; shift 2 ;;
    --lineage) lineage=$2; shift 2 ;;
    --mode) mode=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    --run-name) run_name=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --offline) offline=true; shift ;;
    *) usage; exit 64 ;;
  esac
done

[[ -s "$input" && -n "$lineage" && -n "$output" && -n "$run_name" ]] || { usage; exit 2; }
[[ "$mode" =~ ^(genome|proteins|transcriptome)$ ]] || { echo "Unsupported BUSCO mode: $mode" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "threads must be a positive integer" >&2; exit 2; }
[[ "$run_name" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "run-name contains unsupported characters" >&2; exit 2; }
command -v busco >/dev/null

if [[ -e "$output/$run_name" ]]; then
  echo "Refusing to overwrite an existing BUSCO run: $output/$run_name" >&2
  exit 3
fi
mkdir -p "$output"

busco_args=(
  --in "$input"
  --lineage_dataset "$lineage"
  --mode "$mode"
  --out "$run_name"
  --out_path "$output"
  --cpu "$threads"
)
if [[ "$offline" == true ]]; then
  busco_args+=(--offline)
fi

busco --version > "$output/${run_name}.busco_version.txt"
busco "${busco_args[@]}"
