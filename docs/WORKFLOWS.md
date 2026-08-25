# Workflow guide

## Maintained workflows

| Workflow | Entry point | Configuration | Primary output |
|---|---|---|---|
| Sentieon alignment and gVCF | `python -m eps_workflows.sentieon_gvcf` | `configs/sentieon_gvcf.example.json` | One rendered shell script per sample |
| Sentieon joint SNP calling | `python -m eps_workflows.sentieon_joint_calling` | `configs/sentieon_joint.example.json` plus a gVCF manifest | Joint VCF, SNP VCF and hard-filtered SNP VCF |
| Hi-C scaffolding | `python -m eps_workflows.hic_workflow` | `configs/hic.example.json` | Juicer and 3D-DNA shell script |
| Genome annotation stages | `python -m eps_workflows.annotation_workflow` | `configs/annotation.example.json` | Dependency-ordered stage runner |
| T-DNA analysis | `python pipelines/tdna/run_tdna.py` | `configs/tdna.example.json` | Rendered command or executed pipeline |
| Natural-line mutation rate | `python analysis/mutation_rate/run.py` | `configs/mutation_rate_natural.example.json` plus input manifests | Candidate audit and rate tables |
| Transgenic-line mutation rate | `python analysis/mutation_rate/run.py` | `configs/mutation_rate_transgenic.example.json` plus input manifests | Per-interval candidate and rate tables |
| Simulation benchmark | `python analysis/simulation/benchmark_variants.py` | Command-line inputs | Precision, recall and error sets |
| Sequencing-depth plot | `Rscript analysis/coverage/plot_depth_distribution.R` | Command-line inputs | PDF and PNG plots |
| GO enrichment | `Rscript analysis/enrichment/go_enrichment.R` | Command-line inputs | Enrichment table and GO DAGs |
| Whole-genome Hi-C plot | `python analysis/hic/plot_whole_genome.py` | Command-line inputs | Contact-map figure |
| Final sparse Hi-C heatmap | `Rscript analysis/hic/plot_sparse_contact_heatmap.R` | X/Y/Z matrix, chromosome breaks and labels | PDF and PNG contact heatmaps |
| Five-track SNP Circos | `Rscript analysis/circos/plot_snp_tracks.R` | Five interval tracks and chromosome sizes | PDF, SVG, PNG and TIFF figure |
| Comparative-genome Circos | `Rscript analysis/circos/plot_comparative_genome.R` | TE, gene, LAI, chromosome and synteny tables | PDF, SVG, PNG and TIFF figure |
| Annotation workflow figure | `Rscript analysis/schematics/plot_annotation_workflow.R` | Command-line export settings | PDF, SVG, PNG and TIFF figure |
| SNP workflow figure | `Rscript analysis/schematics/plot_snp_workflow.R` | Window, step and hotspot threshold settings | PDF, SVG, PNG and TIFF figure |
| Conservation remapping | `python analysis/conservation/remap_scores.py` | Command-line inputs | Remapped score table |
| Population VCF table | `python analysis/variants/genotype_table.py` | Multi-sample VCF | Site, genotype and substitution-spectrum tables |
| Hot-region SNP filter | `python analysis/variants/hot_region_filter.py` | VCF and chromosome sizes | Audited windows, merged hot regions and filtered VCF |
| Wide publication variant table | `python analysis/variants/build_publication_table.py` | Filtered VCF and optional ANNOVAR/SnpEff tables | Gzip-ready sample-wide annotated table |
| Variant window density | `python analysis/variants/window_density.py` | Variant table and chromosome sizes | Zero-filled window-density table |
| Chromosome Circos tracks | `python analysis/circos/prepare_tracks.py` | Selected VCF, genes, FASTA and chromosome sizes | Five tracks and karyotype |
| Chromosome Circos renderer | `python analysis/circos/render_config.py` | `configs/circos.example.json` | Sanitized Circos configuration and optional PNG/SVG execution |
| Annotation alignment filter | `python analysis/annotation/filter_blast_by_coverage.py` | Length and BLAST tables | Coverage-filtered alignments |
| SaProt mutation scoring | `python analysis/protein_effect/saprot_mutation_scan.py` | Structure, position and external model installation | Long-form substitution scores |
| Merqury assessment | `workflows/assembly/run_merqury.sh` | Command-line inputs | Merqury quality results |
| NUCmer collinearity | `workflows/comparative/nucmer_collinearity.sh` | Sample manifest | Alignment blocks and chromosome lengths |
| Assembly FASTA statistics | `python analysis/assembly/assembly_stats.py` | Assembly FASTA | Contiguity and composition table |
| Merqury assessment | `workflows/assembly/run_merqury.sh` | One or more read files and an assembly | K-mer completeness and QV results |
| BUSCO assessment | `workflows/assembly/run_busco.sh` | Sequence input, mode and lineage | BUSCO run directory and version record |

## End-to-end pipelines

| Pipeline | Entry point | Composition |
|---|---|---|
| Genome assembly and assessment | `python pipelines/genome_assembly/run.py` | Optional Hi-C rendering, FASTA statistics, Merqury and BUSCO |
| Genome annotation | `python pipelines/genome_annotation/run.py` | Dependency-ordered annotation stages |
| SNP calling | `python pipelines/snp_calling/run.py` | Per-sample Sentieon gVCFs and joint genotyping/hard filtering |
| T-DNA analysis | `python pipelines/tdna/run_tdna.py` | Breakpoint, insertion and copy-number stages |

The genome assembly pipeline starts from a supplied assembly. It does not claim
to reproduce a primary contig-assembler command that was absent from the
audited source roots. BUSCO mode is explicit because the recovered manuscript
job used predicted proteins, whereas BUSCO genome mode answers a different
question.

## Configuration policy

All paths, sample identifiers, tool locations and runtime resources are supplied
through JSON configuration or command-line arguments. Examples use descriptive
placeholders and do not encode a laboratory, user account, host or storage layout.

Configuration files containing real local paths should be stored outside the
repository. Only sanitized examples belong in version control.

## Safety behavior

- Workflow rendering is the default; execution requires an explicit flag.
- The Hi-C renderer does not overwrite regular files or delete existing output
  directories.
- Sentieon intermediate cleanup is opt-in.
- Joint calling uses a manifest rather than embedded sample paths. Its reviewed
  hard-filter expressions are explicit in the JSON configuration; the supplied
  workflow image supports the corrected `HaplotypeScore > 13` default.
- Optional post-dedup BAM filtering is explicit and parameterized. The example
  configuration enables the manuscript MAPQ, SAM-flag, NM and CIGAR rules.
- Annotation stages record completion state and validate declared inputs and
  outputs.
- The T-DNA workflow records stage status and retains risk-labelled candidates
  for review.
- Mutation-rate thresholds and pedigree relationships are explicit in JSON;
  exact Poisson intervals are calculated with SciPy.

## Publication boundary

This repository contains maintained analysis code, synthetic test fixtures and
the curated manuscript figures listed in `figures/SHA256SUMS`. Internal source
scans, historical snapshots, scheduler-generated scripts, software
environments, raw sequencing data and bulk derived data are outside the
public-release boundary.
