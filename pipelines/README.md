# End-to-end pipelines

`pipelines/` contains manuscript-level entry points that compose reusable
workflow stages. Scientific calculations and plot-specific transformations
remain under `analysis/`; reusable external-tool wrappers remain under
`workflows/`; shared Python implementations remain under `src/eps_workflows/`.

Available pipelines:

- `genome_assembly/`: Hi-C script rendering plus assembly statistics, Merqury
  and BUSCO assessment from a supplied assembly;
- `genome_annotation/`: dependency-ordered annotation stages;
- `snp_calling/`: per-sample Sentieon gVCFs followed by joint genotyping; and
- `tdna/`: T-DNA insertion and copy-number analysis.

Rendering is the default. External tools run only when an entry point receives
an explicit `--execute` flag.
