# Hi-C whole-genome contact heatmaps

The final post-review figure source uses a sparse table with `X`, `Y` and `Z`
columns, a numeric chromosome-boundary vector and a chromosome-label file. It
plots log10-transformed contact values, chromosome separators and labels but
does not perform hierarchical clustering or calculate a dendrogram.

Use the parameterized reconstruction for this final matrix format:

```bash
Rscript analysis/hic/plot_sparse_contact_heatmap.R \
  --matrix /path/to/heatmap1000000.txt \
  --breaks /path/to/chromosome_breaks.txt \
  --labels /path/to/chromosome_labels.txt \
  --output-pdf /path/to/HiC_heatmap.pdf \
  --output-png /path/to/HiC_heatmap.png \
  --width 10 \
  --height 10 \
  --dpi 300
```

`analysis/hic/plot_whole_genome.py` remains available for Cooler `.mcool`
inputs. The two scripts represent different matrix formats and should not be
treated as numerically interchangeable without validating bin order,
normalization and resolution.
