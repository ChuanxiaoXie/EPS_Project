# EPS data-processing workflows

This repository contains parameterized, reproducible workflows used for EPS
resequencing and genome analysis. It is prepared for public release with the
associated manuscript.

## Repository contents

- `src/eps_workflows/`: maintained workflow renderers for Sentieon gVCF and
  joint SNP calling, Hi-C scaffolding and genome annotation.
- `pipelines/`: manuscript-level genome assembly/assessment, genome annotation,
  SNP calling and T-DNA entry points.
- `analysis/`: manuscript analyses for assembly statistics, mutation rate, hot-region filtering,
  publication variant tables, Circos, simulation benchmarking, sequencing
  depth, GO enrichment, Hi-C visualization and conservation scores.
- `workflows/`: reusable assembly assessment, annotation, SNP-calling and
  genome-comparison stages.
- `configs/`: example JSON configurations containing placeholders only.
- `docs/WORKFLOWS.md`: workflow scope, inputs, outputs and safety behavior.
- `docs/SOURCE_SELECTION.md`: auditable mapping from legacy script families to
  the public implementations.
- `docs/HOT_REGION_AND_CIRCOS.md`: exact high-density filtering, wide-table and
  five-track Circos analysis definitions.
- `docs/HIC_HEATMAP.md`: distinction between the final sparse-matrix Hi-C
  contact heatmap and the Cooler-based plotting route.
- `docs/FIGURE_WORKFLOWS.md`: R figure inputs, export commands, reviewed
  parameter decisions and visual-QA requirements.
- `scripts/check_public_release.py`: local and CI release-safety scanner.
- `tests/`: unit tests for configuration and workflow rendering.
- `testdata/`: small synthetic fixtures with no real biological observations.
- `figures/`: curated assembly, Hi-C, SNP-density and Circos outputs.
- `docs/PUBLICATION_SCOPE_AUDIT.md`: evidence-based coverage and remaining
  source gaps for the manuscript-facing workflows.

Raw sequencing data, reference genomes, bulk analysis outputs, credentials,
institutional storage paths and internal provenance snapshots are not included.
Only the explicitly curated manuscript figures under `figures/` are retained.

## Quick start

Use Python 3.9 or newer. Rendering workflow scripts requires no third-party
Python packages.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[mutation-rate]"

python -m eps_workflows.sentieon_gvcf \\
  --config configs/sentieon_gvcf.example.json \\
  --script-dir generated/sentieon

python -m eps_workflows.sentieon_joint_calling \\
  --config configs/sentieon_joint.example.json \\
  --gvcf-manifest testdata/snp/gvcfs.manifest.tsv \\
  --script generated/sentieon_joint.sh

python -m eps_workflows.hic_workflow \\
  --config configs/hic.example.json \\
  --script generated/hic.sh

python -m eps_workflows.annotation_workflow \\
  --config configs/annotation.example.json \\
  --script generated/annotation.sh

python pipelines/tdna/run_tdna.py --config configs/tdna.example.json

python pipelines/genome_assembly/run.py \
  --config configs/genome_assembly_pipeline.example.json \
  --script-dir generated/genome_assembly

python pipelines/genome_annotation/run.py \
  --config configs/annotation.example.json \
  --script generated/genome_annotation.sh

python pipelines/snp_calling/run.py \
  --gvcf-config configs/sentieon_gvcf.example.json \
  --joint-config configs/sentieon_joint.example.json \
  --output-dir generated/snp_calling

python analysis/mutation_rate/run.py \
  --config configs/mutation_rate_natural.example.json \
  --gvcf-manifest configs/mutation_rate_natural.gvcfs.example.tsv \
  --bam-manifest configs/mutation_rate_natural.bams.example.tsv \
  --joint-vcf /path/to/joint.vcf.gz \
  --reference /path/to/reference.fa \
  --mask-bed /path/to/error_mask.bed \
  --output /path/to/mutation_rate_output
```

Rendering is the default. Add `--execute` only after replacing every
placeholder and reviewing the generated shell script.

## Validation

The R figure smoke tests require `ggplot2`, `RColorBrewer`, `circlize` and
`svglite`.

```bash
python scripts/check_testdata_coverage.py
python scripts/run_example_smoke_tests.py
python scripts/check_figure_hashes.py
python -m unittest discover -s tests -v
python scripts/check_public_release.py --root .
python -m compileall -q src pipelines analysis scripts
bash -n pipelines/tdna/pipeline.sh
bash -n workflows/assembly/run_merqury.sh
bash -n workflows/assembly/run_busco.sh
bash -n workflows/comparative/nucmer_collinearity.sh
Rscript -e 'for (f in list.files("analysis", pattern="[.]R$", recursive=TRUE, full.names=TRUE)) parse(f)'
```

The test-data coverage check audits the fixture mapping for every public
analysis and pipeline entry point. The smoke test then executes lightweight
analyses, renders workflows that depend on external domain tools, and creates
the available R figures. See `testdata/test_matrix.tsv` for each entry point's
fixture paths and exact test mode.

The release scanner rejects private storage roots, personal home directories,
network identifiers, email addresses, credential-like strings, private-key
material, sensitive filenames and raw-data files. Important code comments are
written in English.
