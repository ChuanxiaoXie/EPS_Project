#!/usr/bin/env bash
set -Eeuo pipefail

if (($# != 4)); then
  echo "Usage: $0 REFERENCE_FASTA SAMPLES_TSV OUTPUT_DIR THREADS" >&2
  echo "SAMPLES_TSV columns: sample_id<TAB>assembly_fasta" >&2
  exit 64
fi

reference=$1
samples_tsv=$2
output_root=$3
threads=$4
[[ -s "$reference" && -s "$samples_tsv" && "$threads" =~ ^[1-9][0-9]*$ ]] || exit 2
for tool in nucmer delta-filter show-coords samtools; do command -v "$tool" >/dev/null; done
mkdir -p "$output_root/reference"
samtools faidx "$reference"
awk 'BEGIN{OFS="\t"}{print "reference",$1,$2}' "$reference.fai" > "$output_root/reference/chromosome_lengths.tsv"

while IFS=$'\t' read -r sample assembly extra; do
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  [[ -z "${extra:-}" && "$sample" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ && -s "$assembly" ]] || {
    echo "Invalid sample row: $sample" >&2
    exit 2
  }
  sample_dir="$output_root/$sample"
  mkdir -p "$sample_dir"
  nucmer -t "$threads" --mum -p "$sample_dir/$sample" "$reference" "$assembly"
  delta-filter -i 90 -l 15000 -q "$sample_dir/$sample.delta" > "$sample_dir/$sample.filtered.delta"
  show-coords -THrd "$sample_dir/$sample.filtered.delta" > "$sample_dir/$sample.coords.tsv"
  awk -v sample="$sample" 'BEGIN {
      OFS="\t";
      print "block_id","reference","reference_chrom","reference_start","reference_end","query","query_chrom","query_start","query_end","strand"
    }
    {
      ref_start=($1<$2)?$1:$2; ref_end=($1>$2)?$1:$2;
      query_start=($3<$4)?$3:$4; query_end=($3>$4)?$3:$4;
      strand=($9=="-1")?"-":"+";
      print NR,"reference",$10,ref_start,ref_end,sample,$11,query_start,query_end,strand
    }' "$sample_dir/$sample.coords.tsv" > "$sample_dir/$sample.blocks.tsv"
  samtools faidx "$assembly"
  awk -v sample="$sample" 'BEGIN{OFS="\t"}{print sample,$1,$2}' "$assembly.fai" > "$sample_dir/$sample.chromosome_lengths.tsv"
done < "$samples_tsv"
