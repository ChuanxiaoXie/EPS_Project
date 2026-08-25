# SNP-calling workflows

Reusable per-sample gVCF and joint-calling renderers are implemented in:

- `src/eps_workflows/sentieon_gvcf.py`; and
- `src/eps_workflows/sentieon_joint_calling.py`.

The end-to-end composition lives in `pipelines/snp_calling/run.py`. Downstream
hot-region filtering and figure preparation remain in `analysis/` so that
variant generation and manuscript analysis are not conflated.
