# Publication-scope coverage audit

This audit covers the two supplied project roots and the current public tree.
It distinguishes a reproducible script from an output-only figure and does not
claim that a workflow exists when the source evidence is absent.

| Area | Public scripts | Curated figure | Status |
|---|---|---|---|
| Primary contig assembly | None located | None | **Source gap.** Across 669 genome-side code files, no Hifiasm, Flye, Canu, NextDenovo, FALCON, purge_dups or equivalent primary-assembly invocation was found. The supplied genome root starts with Hi-C/scaffolding activities. |
| Hi-C chromosome scaffolding | `src/eps_workflows/hic_workflow.py` | `figures/hic/EPSV2_juicerbox.pdf` | Covered for Juicer and 3D-DNA script rendering; external tools and real reads are not bundled. |
| Final Hi-C contact heatmap | `analysis/hic/plot_sparse_contact_heatmap.R` | `figures/hic/HiC_heatmap_adjusted_assembly.pdf` | Covered. The source uses sparse X/Y/Z contacts and chromosome break/label files. It is not a hierarchical clustering heatmap. |
| Assembly and annotation quality | `analysis/assembly/assembly_stats.py`, `workflows/assembly/run_merqury.sh`, `workflows/assembly/run_busco.sh` | `figures/assembly/busco_figure.pdf` | Merqury paired short-read and single long-read jobs were recovered and unified. The recovered BUSCO job used predicted proteins, embryophyta_odb10, 15 threads and offline mode; its software version still requires confirmation. |
| Per-sample SNP calling | `src/eps_workflows/sentieon_gvcf.py` | `figures/reference/author_reference_snp_workflow.png` | Covered by parameterized Sentieon alignment, duplicate removal, reviewed post-dedup BAM filtering and Haplotyper gVCF rendering. |
| Joint SNP calling and hard filtering | `src/eps_workflows/sentieon_joint_calling.py` | `figures/reference/author_reference_snp_workflow.png` | Covered by Sentieon GVCFtyper, SNP selection and the reviewed hard-filter expressions. |
| High-density-region filtering and publication table | `analysis/variants/hot_region_filter.py`, `build_publication_table.py`, `genotype_table.py` | `figures/snp_density/` | Covered, including the 200 kb window, 20 kb step and count >100 rule used by the final source. |
| Five-track chromosome Circos | `analysis/circos/prepare_tracks.py`, `render_config.py`, `plot_snp_tracks.R` | `figures/circos/circos_200k_0727_fixTH_scatter_abcde.png`, `figures/reference/author_reference_snp_circos_panel.png` | Covered for track preparation, direct a-e track labels, original Circos configuration and the new R renderer. |
| Comparative-genome Circos | `analysis/circos/plot_comparative_genome.R`, `workflows/comparative/nucmer_collinearity.sh` | `figures/reference/author_reference_supplementary_figure4.png` | The missing plotting layer was reconstructed in R from TE, gene, LAI and collinearity input tables. |
| Annotation workflow diagram | `analysis/schematics/plot_annotation_workflow.R` | `figures/reference/author_reference_supplementary_figure4.png` | Covered as a publication schematic; quantitative annotation commands remain in the dependency-ordered workflow renderer. |

## Test boundary

`testdata/` contains a synthetic two-chromosome assembly, two-sample VCF,
annotation-shaped tables and deterministic ten-chromosome R figure inputs.
Running `python scripts/run_example_smoke_tests.py` executes dependency-free
analysis steps, validates the three top-level biological pipelines and three
external-workflow renderers, and renders PDF and PNG outputs for all six R
figure workflows when their packages are available.

The synthetic fixtures prove that the maintained code paths execute; they do
not replace external-tool integration tests for Sentieon, GATK, Juicer,
3D-DNA, Circos, Merqury or BUSCO.

## Author confirmations still required

1. Supply the primary contig-assembly command, software version and key
   parameters, or confirm that assembly was delivered externally and should be
   described without a public launch script.
2. Confirm the BUSCO software version and that the recovered protein-mode
   assessment is the result intended for the included completeness figure.
3. Confirm in the final Methods text that the reviewed public default
   `HaplotypeScore > 13.0` is intended. The retained legacy script encoded the
   opposite operator, while the author-supplied workflow image encodes `>`.
