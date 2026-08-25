# Genome annotation workflow

The reusable annotation-stage renderer is implemented in
`src/eps_workflows/annotation_workflow.py` and configured by
`configs/annotation.example.json`. The manuscript-level entry point is
`pipelines/genome_annotation/run.py`.

This separation keeps one maintained implementation while making both the
stage workflow and the end-to-end pipeline easy to discover.
