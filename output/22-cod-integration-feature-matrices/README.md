# output/22-cod-integration-feature-matrices

Outputs from [`code/22-cod-integration-feature-matrices.Rmd`](../../code/22-cod-integration-feature-matrices.Rmd).

This notebook prepares **matched feature matrices** for the paired samples:
normalized gene expression (from the `06.2` featureCounts gene matrix, normalized with
DESeq2 VST / edgeR logCPM) and gene-level (gene-body and promoter) DNA methylation
summaries derived from the methylseq `coverage2cytosine` files.

Expected files (created when the notebook is run):

- `expression-vst-paired.csv` — VST-normalized gene expression, paired samples, genes x samples.
- `methylation-genebody-paired.csv` — gene-body weighted methylation, paired samples.
- `methylation-promoter-paired.csv` — promoter (TSS-flank) weighted methylation, paired samples.
- `common-genes.txt` — gene IDs present in both expression and methylation matrices.
- `feature-matrix-summary.csv` — per-matrix dimensions and coverage/filtering notes.

This is a scaffold/draft. Run the notebook to populate this directory.
