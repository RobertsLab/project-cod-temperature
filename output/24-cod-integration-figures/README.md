# output/24-cod-integration-figures

Outputs from [`code/24-cod-integration-figures.Rmd`](../../code/24-cod-integration-figures.Rmd).

This notebook renders **publication-quality figures** for the integrated RNA-seq + WGBS
analysis, saved at high DPI (300+).

Expected files (created when the notebook is run):

- `fig-joint-pca.png` / `.pdf` — joint multi-omics PCA, colored by `temp_treatment`.
- `fig-expr-meth-volcano.png` / `.pdf` — per-gene correlation vs -log10(FDR).
- `fig-top-gene-scatter.png` / `.pdf` — expression-vs-methylation scatter for top genes (faceted).
- `fig-integration-heatmap.png` / `.pdf` — heatmap of top integrated features across paired samples.

This is a scaffold/draft. Run the notebook to populate this directory.
