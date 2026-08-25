# Synthetic mutation-rate fixtures

The opportunity and candidate tables directly exercise final mutation-rate
calculation when SciPy is installed. The configuration and manifests also
render the complete command plan.

The small gVCF-shaped text files document the expected record structure. Full
callable-genome and candidate-discovery integration tests require indexed
gVCFs, coordinate-sorted BAMs, bcftools, bedtools and samtools. Binary files
with false `.bam` or `.g.vcf.gz` labels are intentionally not included.
