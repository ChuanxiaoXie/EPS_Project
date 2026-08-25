# SNP calling pipeline

This pipeline composes per-sample Sentieon alignment/gVCF scripts and the joint
GVCFtyper plus reviewed SNP hard-filter script. It writes a generated gVCF
manifest and one master script without executing external software by default.

```bash
python pipelines/snp_calling/run.py \
  --gvcf-config configs/sentieon_gvcf.example.json \
  --joint-config configs/sentieon_joint.example.json \
  --output-dir generated/snp_calling
```

Downstream high-density-region filtering, publication tables and Circos tracks
remain under `analysis/variants/` and `analysis/circos/` because they consume
the called variants rather than producing the primary joint call set.
