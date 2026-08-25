#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 --read FILE [--read FILE ...] --assembly FILE --output DIR --merqury-root DIR [--prefix NAME] [--kmer INT]" >&2
}

reads=()
assembly=""
output=""
merqury_root=""
prefix="assembly"
kmer=""
while (($#)); do
  case "$1" in
    --read) reads+=("$2"); shift 2 ;;
    --read1) reads+=("$2"); shift 2 ;;
    --read2) reads+=("$2"); shift 2 ;;
    --assembly) assembly=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    --merqury-root) merqury_root=$2; shift 2 ;;
    --prefix) prefix=$2; shift 2 ;;
    --kmer) kmer=$2; shift 2 ;;
    *) usage; exit 64 ;;
  esac
done

[[ ${#reads[@]} -gt 0 && -s "$assembly" && -d "$merqury_root" && -n "$output" ]] || {
  usage
  exit 2
}
for read_file in "${reads[@]}"; do
  [[ -s "$read_file" ]] || { echo "Missing read file: $read_file" >&2; exit 2; }
done
command -v seqkit >/dev/null
command -v meryl >/dev/null
[[ -x "$merqury_root/best_k.sh" && -x "$merqury_root/merqury.sh" ]] || {
  echo "Invalid Merqury installation: $merqury_root" >&2
  exit 2
}
if [[ -d "$output" && -n "$(find "$output" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to reuse a non-empty output directory: $output" >&2
  exit 3
fi
mkdir -p "$output"
cd "$output"

assembly_length=$(seqkit stat -T "$assembly" | awk -F '\t' 'NR==2 {print $5}')
k=${kmer:-$($merqury_root/best_k.sh "$assembly_length" | awk 'NR==3 {print int($0 + 0.5)}')}
[[ "$k" =~ ^[0-9]+$ ]] || { echo "Unable to determine k-mer size" >&2; exit 4; }

read_databases=()
for index in "${!reads[@]}"; do
  database="read_$((index + 1)).meryl"
  meryl k="$k" count output "$database" "${reads[$index]}"
  read_databases+=("$database")
done
if [[ ${#read_databases[@]} -eq 1 ]]; then
  meryl union-sum output reads.meryl "${read_databases[0]}"
else
  meryl union-sum output reads.meryl "${read_databases[@]}"
fi
"$merqury_root/merqury.sh" reads.meryl "$assembly" "$prefix"
