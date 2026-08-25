# Genome assembly and assessment pipeline

This pipeline starts from a supplied assembly because no primary contig-
assembler command was recovered. It can render an optional Hi-C scaffolding
script and then run FASTA statistics, Merqury and BUSCO against explicitly
configured assessment inputs.

```bash
python pipelines/genome_assembly/run.py \
  --config configs/genome_assembly_pipeline.example.json \
  --script-dir generated/genome_assembly
```

Review every generated script before adding `--execute`. If Hi-C is enabled,
`assessment_assembly` must point to the intended post-scaffolding FASTA output.
BUSCO protein mode evaluates the predicted gene set; genome mode evaluates the
assembly sequence. These claims must not be interchanged in the manuscript.
