`project-cod-temperature/code`

Directory containing scripts, primarily in R Markdown. Numbered files have a corresponding output directory in `../output`.

R Markdown files may have been knitted to HTML and/or markdown. These files have suffixes of `.html` or `.md`.

---

- [`13.0.0-RNAseq-edgeR.Rmd`](./13.0.0-RNAseq-edgeR.Rmd): edgeR RNA-seq analysis to identify differentially expressed genes between 9<sup>o</sup>C and 16<sup>o</sup>C samples. Output files are in [`../output/13.0.0-RNAseq-edgeR/`](../output/13.0.0-RNAseq-edgeR/).

- [`14.00-heatwave-genetics-trimmed-bwa-alignments.Rmd`](./14.00-heatwave-genetics-trimmed-bwa-alignments.Rmd): BWA alignment of trimmed reads from the heatwave genetics project. Output files are in [`../output/14.00-heatwave-genetics-trimmed-bwa-alignments/`](../output/14.00-heatwave-genetics-trimmed-bwa-alignments/).

---

### Integrated RNA-seq + DNA methylation series (`21`–`24`)

A scaffolded, end-to-end workflow that joins the RNA-seq and WGBS methylation data for the paired-sample subset (individuals assayed by both platforms) and carries it from raw inputs to publication figures. Run in order; each notebook consumes the previous one's outputs.

- [`21-cod-integration-crosswalk.Rmd`](./21-cod-integration-crosswalk.Rmd): Builds the sample crosswalk / integrated metadata table linking RNA-seq `sample_{N}`, WGBS `{NN}B`, and `GeneticSamplingCount`, and flags the paired-sample subset. Output files are in [`../output/21-cod-integration-crosswalk/`](../output/21-cod-integration-crosswalk/).

- [`22-cod-integration-feature-matrices.Rmd`](./22-cod-integration-feature-matrices.Rmd): Prepares matched feature matrices for paired samples — DESeq2/edgeR-normalized gene expression and gene-body/promoter methylation summarized from the methylseq coverage files. Output files are in [`../output/22-cod-integration-feature-matrices/`](../output/22-cod-integration-feature-matrices/).

- [`23-cod-integration-analysis.Rmd`](./23-cod-integration-analysis.Rmd): Expression–methylation integration — per-gene correlation (with FDR), joint PCA, and optional mixOmics DIABLO multi-omics model. Output files are in [`../output/23-cod-integration-analysis/`](../output/23-cod-integration-analysis/).

- [`24-cod-integration-figures.Rmd`](./24-cod-integration-figures.Rmd): Publication-quality figures (joint PCA, correlation volcano, top-gene scatters, integration heatmaps) saved at high DPI. Output files are in [`../output/24-cod-integration-figures/`](../output/24-cod-integration-figures/).
