# Synthetic R-figure fixtures

`generate_r_figure_fixtures.R` creates deterministic ten-chromosome tracks for
the five-track SNP Circos plot, the CA7301/B73 comparative Circos plot and the
whole-genome Hi-C contact heatmap. The generated values are artificial and are
only suitable for executable rendering tests.

```bash
Rscript testdata/figures/generate_r_figure_fixtures.R /tmp/eps-figure-fixtures
```
