# Source selection and refactoring record

## Selection basis

The internal discovery scan inspected 121,253 code files and classified 5,422
as project-authored source candidates. The public repository does not reproduce
that raw tree. It retains one maintained implementation for each manuscript
method while excluding generated scripts, scheduler copies, sample-expanded
duplicates, report assets, software environments and temporary files.

Selection favored scripts that encode a scientific threshold, formula,
candidate definition, coordinate transformation or final figure. Literal job
submission wrappers and files that only substitute sample paths were replaced
by configuration-driven implementations.

## Public mapping

| Legacy analysis family | Public implementation | Refactoring decision |
|---|---|---|
| Repeated Sentieon FASTQ-to-gVCF jobs | `src/eps_workflows/sentieon_gvcf.py` | Consolidated into one per-sample renderer. |
| Sentieon GVCFtyper jobs plus final SNP hard filtering | `src/eps_workflows/sentieon_joint_calling.py` | Replaced sample-expanded `-v` arguments with a manifest and made every hard-filter expression explicit. |
| Genome Hi-C, Juicer and 3D-DNA jobs | `src/eps_workflows/hic_workflow.py` | Parameterized and made non-destructive. |
| Genome-annotation stage wrappers | `src/eps_workflows/annotation_workflow.py` | Replaced by a dependency-aware runner. |
| T-DNA assembly and breakpoint scripts | `pipelines/tdna/` | Retained core algorithms and parameterized all runtime inputs. |
| Natural and transgenic lineage mutation-rate scripts | `analysis/mutation_rate/` | Preserved callable-region rules, first-appearance logic, BAM validation, lineage checks and Poisson intervals; pedigree structure moved to JSON. |
| Sequencing-depth figure copies | `analysis/coverage/plot_depth_distribution.R` | Replaced repeated path-specific R files with one argument-driven plotter. |
| Simulated-read variant comparisons | `analysis/simulation/benchmark_variants.py` | Replaced repeated depth/sample wrappers with a normalized VCF truth benchmark. |
| GOseq and topGO scripts | `analysis/enrichment/go_enrichment.R` | Replaced working-directory state and sourced local objects with explicit inputs. |
| Whole-genome Hi-C plotting | `analysis/hic/plot_whole_genome.py` | Added command-line resolution, chromosome and color-scale parameters. |
| `07.post-review/Readjust/heatmap.juice.R` | `analysis/hic/plot_sparse_contact_heatmap.R` | Reconstructed the final X/Y/Z contact heatmap with explicit matrix, chromosome-boundary, label and output arguments; removed private paths and non-English comments. |
| Conservation-coordinate remapping | `analysis/conservation/remap_scores.py` | Removed fixed filenames and made missing-score behavior explicit. |
| Population VCF genotype statistics | `analysis/variants/genotype_table.py` | Replaced fixed FORMAT-column offsets with FORMAT-key parsing and long-form output. |
| `24.V2Genome.AddM4GAAB/HotRegion.Original` high-density filter | `analysis/variants/hot_region_filter.py` | Preserved the 200-kb/20-kb sliding windows, strict count >100 rule, merged regions and PASS-inside/all-outside SNP retention logic. |
| `popSNP_GQDP_etal.Stat.index.v3_425.pl` plus ANNOVAR/SnpEff join | `analysis/variants/build_publication_table.py` | Rebuilt the actual sample-wide table, removed fixed FORMAT offsets and made the legacy inner join optional. |
| Variant-density window scripts | `analysis/variants/window_density.py` | Consolidated repeated chromosome loops into a zero-filled, group-aware implementation. |
| `HotRegion.Original/Circos/circos_0727_200k` | `analysis/circos/` and `configs/circos.example.json` | Reconstructed all five final 200-kb tracks and the later display-range preset; removed symlinked private inputs and external environment paths. |
| Transcript-support BLAST filtering | `analysis/annotation/filter_blast_by_coverage.py` | Preserved the 80% identity/coverage defaults and exposed both thresholds. |
| SaProt single-position mutation scoring | `analysis/protein_effect/saprot_mutation_scan.py` | Replaced model installation paths, protein naming and device selection with arguments. |
| Paired short-read and single long-read Merqury assessment | `workflows/assembly/run_merqury.sh` | Unified both recovered jobs behind repeated `--read` arguments; removed user environments, email metadata and destructive rerun behavior. |
| BUSCO protein completeness assessment | `workflows/assembly/run_busco.sh` | Recovered the confirmed proteins/embryophyta_odb10/offline/15-thread definition while parameterizing mode, lineage and paths. |
| NUCmer collinearity jobs | `workflows/comparative/nucmer_collinearity.sh` | Replaced personal Perl/tool paths with standard executables and a sample manifest. |

## Scientific invariants retained

- Default callable bases require depth 5–40 and genotype quality at least 30.
- Strict SNV sites use QD at least 2, FS at most 40, mapping quality at least
  50, MQRankSum at least -12.5 and ReadPosRankSum at least -4.
- Target alternate-read support is at least five reads; heterozygous allele
  balance is 0.25–0.75 and homozygous-alternate balance is at least 0.90.
- Strict BAM confirmation requires forward and reverse alternate-read support.
- Mutation-rate denominators are callable base pairs multiplied by generations;
  confidence intervals are exact Poisson intervals.
- Hot regions use 200-kb windows at 20-kb steps and a strict SNP-count rule of
  greater than 100; only PASS SNPs are retained inside merged hot regions.
- The retained legacy source encoded `haplotype_score < 13.0`, while the
  author-supplied workflow image encodes the upper-tail filter. The public
  renderer defaults to `haplotype_score > 13.0`, exposes the operator in JSON
  and records the discrepancy for final Methods confirmation.
- The final Circos data tracks use 200-kb non-overlapping windows for all SNPs,
  genic SNPs, gene features, G+C bases and G/C-to-A/T SNPs.

All values are visible and overridable in the public JSON configurations. A
change to these defaults changes the analysis definition and should be reported
as a sensitivity analysis.

## Excluded after review

- Bundled report JavaScript, software environments and model source trees are
  third-party dependencies, not manuscript analysis code.
- Sample-expanded jobs, scheduler copies and dated reruns are represented by
  their configuration-driven workflow renderers.
- Threshold experiments and incomplete prototypes were replaced by maintained
  implementations with explicit defaults and validation.
- Primer-design and other experimental-support utilities were not part of the
  upstream biological data-processing workflow and remain outside this package.
