# Curated figure outputs

This directory contains only manuscript-relevant derived figures that were
matched to a maintained public analysis or to a documented assessment result.
No raw sequencing data, sample-level source records or private paths are
embedded in the filenames.

## Included groups

- `assembly/`: BUSCO completeness figure supplied as PDF plus a PNG preview.
  The recovered launch script used predicted proteins, embryophyta_odb10,
  15 threads and offline mode. BUSCO software version and exact figure/run
  correspondence still require confirmation from the authors' Methods records.
- `hic/`: the final adjusted-assembly chromosome-partitioned contact heatmap
  and the EPSV2 supplementary Juicer-style contact map. These are contact
  heatmaps, not hierarchical clustering plots, and contain no dendrogram.
- `circos/`: the final 200 kb `fixTH` Circos image with direct a-e track labels,
  matched to the five-track public preparation and configuration scripts.
- `snp_density/`: six original SNP-count distribution plots used to compare
  window/step choices, including the 200 kb/20 kb setting used by the hot-region
  analysis.
- `reference/`: three author-supplied preview images used to reconstruct and
  visually review the R plotting scripts. They are reference screenshots, not
  final high-resolution submission files or synthetic test outputs.

The original outputs and author reference previews are not regenerated during
the small synthetic smoke test. The smoke test writes temporary figures only.

`SHA256SUMS` records the exact release copies. Validate them with
`python scripts/check_figure_hashes.py`.
