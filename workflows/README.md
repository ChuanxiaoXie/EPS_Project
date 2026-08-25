# Reusable workflows

`workflows/` contains reusable stage-level operations used by one or more
top-level pipelines. It is intentionally not a second copy of the Python
implementation under `src/eps_workflows/`.

- `assembly/`: Merqury and BUSCO assessment wrappers;
- `annotation/`: annotation-stage workflow contract and entry points;
- `snp_calling/`: Sentieon gVCF and joint-calling workflow contract; and
- `comparative/`: NUCmer-based collinearity.
