# Genome annotation pipeline

This entry point exposes the dependency-ordered annotation pipeline at the
manuscript level while keeping its reusable implementation in
`src/eps_workflows/annotation_workflow.py`.

```bash
python pipelines/genome_annotation/run.py \
  --config configs/annotation.example.json \
  --script generated/genome_annotation.sh
```

The example stages are placeholders. A publication run must explicitly define
repeat annotation, structural annotation, non-coding RNA annotation,
transcript-supported model training, evidence integration and functional
annotation commands with declared inputs and outputs.
