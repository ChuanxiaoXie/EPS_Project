# Hot-region filtering, publication table and Circos plot

## Provenance

These public implementations reconstruct the final analysis family under
`13.FilterBam/24.V2Genome.AddM4GAAB/HotRegion.Original`. The source family
contains the high-density-region filter, the 78-sample wide variant table and
the final `circos_0727_200k` five-track chromosome plot.

The public code replaces fixed storage paths, dated output names, sample-count
assumptions and positional FORMAT parsing. It does not include source VCFs,
reference sequences, annotation tables or generated tracks. Curated final and
window-comparison figures are retained under `figures/` with checksums.

## Hot-region retention rule

The manuscript workflow counts SNPs in 200-kb windows at 20-kb steps. Windows
with more than 100 SNPs are merged. SNPs inside the merged high-density regions
must have `FILTER=PASS`; SNPs outside those regions are retained regardless of
their FILTER value. Only SNP records are written to the filtered VCF.

```bash
python analysis/variants/hot_region_filter.py \
  --vcf /path/to/merged.vcf.gz \
  --chrom-sizes /path/to/chrom.sizes \
  --output-vcf /path/to/results/hot_region_filtered.vcf.gz \
  --hot-regions /path/to/results/hot_regions.bed \
  --window-table /path/to/results/sliding_window_counts.tsv \
  --positions-output /path/to/results/retained_positions.tsv
```

The window table records every count and decision, making the threshold fully
auditable. `--window-size`, `--step-size`, `--hot-count-threshold` and repeated
`--accepted-hot-filter` arguments expose the analysis definition.

## Wide publication table

The legacy workflow expanded each sample into genotype, total depth, reference
depth, alternate depth, genotype quality and alternate-depth fraction columns.
It then added ANNOVAR regional/gene/protein-change fields and SnpEff effect and
amino-acid fields. The maintained implementation reads FORMAT keys by name and
supports phased, multi-allelic and reordered FORMAT values.

```bash
python analysis/variants/build_publication_table.py \
  --vcf /path/to/results/hot_region_filtered.vcf.gz \
  --annovar-variant-function /path/to/variant_function \
  --annovar-exonic-variant-function /path/to/exonic_variant_function \
  --snpeff-tsv /path/to/snpeff_annotations.tsv \
  --require-annovar \
  --output /path/to/results/publication_variant_table.tsv.gz
```

`--require-annovar` reproduces the legacy inner-join behavior. Omitting it is
safer for new analyses because unannotated variants remain in the table with
explicit `NA` fields.

## Five Circos tracks

The final plot uses non-overlapping 200-kb windows and five data tracks:

1. all selected SNPs;
2. selected SNPs overlapping gene intervals;
3. gene/transcript feature density;
4. reference G+C base count; and
5. G/C-to-A/T SNP density.

Prepare the tracks and karyotype:

```bash
python analysis/circos/prepare_tracks.py \
  --vcf /path/to/selected_snps.vcf.gz \
  --genes-bed /path/to/gene_intervals.bed \
  --reference /path/to/reference.fa \
  --chrom-sizes /path/to/chrom.sizes \
  --window-size 200000 \
  --output-dir /path/to/circos_tracks
```

Render the sanitized Circos configuration, review it, and then execute:

```bash
python analysis/circos/render_config.py \
  --config configs/circos.example.json \
  --output-config /path/to/circos_tracks/circos.conf

circos -conf /path/to/circos_tracks/circos.conf \
  -png -svg \
  -outputdir /path/to/circos_figure \
  -outputfile eps_chromosome_circos
```

The example configuration records the later display-range preset used by the
final source version: 0–65, 0–25, 0–15, 0–119,408 and 0–49 for the five tracks.
These are visualization limits, not variant-selection thresholds. The earlier
source preset used 0–80, 0–43, 0–299, 7,507–119,408 and 0–49.

## Coordinate compatibility

The legacy shell scripts wrote VCF `POS` directly as a BED start. The public
default correctly converts one-based VCF positions to zero-based BED positions.
For a byte-level comparison with legacy boundary behavior, pass
`--legacy-pos-as-bed-start` to the hot-region or Circos track command and report
that compatibility mode in the analysis record.
