# Mutation-rate workflow

## Inputs

- A JSON pedigree configuration.
- A two-column `sample`/`gvcf` manifest.
- A two-column `sample`/`bam` manifest in the same biological sample order.
- A joint VCF containing all configured samples.
- The indexed reference FASTA and an error-mask BED.

## Stages

1. Convert each gVCF into merged callable regions using the configured depth and
   genotype-quality thresholds.
2. Intersect the samples required for each pedigree interval and subtract the
   error mask.
3. Identify biallelic first-appearance SNVs and enforce ancestral,
   transmission, site-quality and allele-balance rules.
4. Re-query candidate loci from the ordered BAM files with mapping quality 50
   and base quality 30 by default.
5. Validate strand support and lineage consistency from mpileup evidence.
6. Exclude sites that first appear independently in more than one configured
   lineage from the primary independent-event numerator.
7. Divide accepted independent events by callable site-generations and report
   exact Poisson 95% confidence intervals.

The natural and transgenic examples each include matching gVCF and BAM
manifests under `configs/`. Replace only the placeholder paths; sample names
must remain consistent with the selected pedigree JSON.

`analysis/mutation_rate/run.py` prints the complete command plan by default.
Add `--execute` only after checking every input and sample mapping. The two
example configurations preserve the natural-line and transgenic-line pedigree
structures used by the analysis while keeping all storage locations external.
