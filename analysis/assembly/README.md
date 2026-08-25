# Assembly analysis

This module covers the assembly analyses that can be reconstructed from the
audited source material:

- streaming FASTA contiguity and composition statistics;
- Merqury k-mer completeness and consensus-quality assessment; and
- BUSCO completeness assessment and the retained BUSCO figure.

`plot_busco_summary.R` reads BUSCO category counts from a TSV source table;
the manuscript counts are stored beside the curated figure rather than being
hard-coded in the plotting script.

The original project contains paired short-read and single long-read Merqury
jobs. The maintained `workflows/assembly/run_merqury.sh` accepts either one or
multiple `--read` arguments and does not delete previous results.

The recovered BUSCO job evaluated predicted proteins with
`embryophyta_odb10`, 15 CPUs and offline mode. It therefore measures annotation
completeness rather than directly proving primary contig assembly quality. The
public BUSCO workflow exposes the mode and lineage explicitly.

No primary contig-assembler launch command was found in the supplied project
roots. The public genome pipeline consequently starts from an existing draft
assembly and covers Hi-C scaffolding plus assembly assessment. It must not be
cited as evidence for an unobserved Hifiasm, Flye, Canu or similar run.

```bash
python analysis/assembly/assembly_stats.py \
  --assembly /path/to/assembly.fa \
  --output /path/to/assembly_stats.tsv

bash workflows/assembly/run_merqury.sh \
  --read /path/to/read1.fastq.gz \
  --read /path/to/read2.fastq.gz \
  --assembly /path/to/assembly.fa \
  --output /path/to/merqury \
  --merqury-root /path/to/merqury

bash workflows/assembly/run_busco.sh \
  --input /path/to/predicted_proteins.fa \
  --lineage /path/to/embryophyta_odb10 \
  --mode proteins \
  --output /path/to/busco \
  --run-name assembly_busco \
  --threads 15 \
  --offline

Rscript analysis/assembly/plot_busco_summary.R \
  --summary figures/assembly/busco_summary.tsv \
  --output-prefix /path/to/busco_summary
```
