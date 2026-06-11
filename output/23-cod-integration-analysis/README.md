# output/23-cod-integration-analysis

Outputs from [`code/23-cod-integration-analysis.Rmd`](../../code/23-cod-integration-analysis.Rmd).

This notebook runs the **expression–methylation integration analysis** on the paired
samples: per-gene expression-vs-methylation correlation (with FDR), a joint PCA of the
combined omics blocks, and an optional multi-omics integration (mixOmics/DIABLO).

Expected files (created when the notebook is run):

- `gene-expr-meth-correlation.csv` — per-gene correlation coefficient, p-value, BH FDR.
- `joint-pca-scores.csv` — sample scores from PCA on the concatenated/scaled blocks.
- `joint-pca-variance.csv` — variance explained per component.
- `diablo-loadings.csv` — (optional) mixOmics DIABLO feature loadings, if run.
- `integration-session-info.txt` — `sessionInfo()` for reproducibility.

This is a scaffold/draft. Run the notebook to populate this directory.
