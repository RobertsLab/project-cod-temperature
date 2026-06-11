# output/21-cod-integration-crosswalk

Outputs from [`code/21-cod-integration-crosswalk.Rmd`](../../code/21-cod-integration-crosswalk.Rmd).

This notebook builds the **sample crosswalk / integrated metadata table** that links
RNA-seq sample IDs (`sample_{GeneticSamplingCount}`), WGBS methylation IDs (`{NN}B`),
and the shared `GeneticSamplingCount` key, and flags the paired-sample subset that
carries both assays.

Expected files (created when the notebook is run):

- `integrated-sample-metadata.csv` — full crosswalk (all RNA-seq and methylation samples).
- `paired-sample-metadata.csv` — paired subset only (both assays present).
- `crosswalk-summary.csv` — counts by assay availability and `temp_treatment`.

This is a scaffold/draft. Run the notebook to populate this directory.
