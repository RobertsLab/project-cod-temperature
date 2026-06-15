# Integrated Methylation and RNA-seq analysis

Goal: integrate blood RNA-seq and DNA methylation (WGBS) from juvenile Pacific cod held
at two temperatures to understand **epigenetic control of gene expression** under thermal
stress.

## Study design

- **Genome:** NCBI RefSeq `GCF_031168955.1`
  (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_031168955.1/)
- **Warm group:** samples 1–36 = blood from cod held at **16 °C**
- **Cold group:** samples 44–74 = blood from cod held at **0 °C**
- **Key:** WGBS and RNA-seq use the **same `{NN}B` sample IDs for the same individuals**
  (e.g. `13B` is one fish assayed by both), so integration is per-individual.

## Raw data (Owl)

- WGBS: http://owl.fish.washington.edu/nightingales/G_macrocephalus/AN00025267/
- RNA-seq: http://owl.fish.washington.edu/nightingales/G_macrocephalus/AN00025268/

## Analysis plan (notebooks 25–31)

| # | Notebook | Purpose |
|---|----------|---------|
| 25 | [`25-blood-sample-sheet.Rmd`](25-blood-sample-sheet.Rmd) | Discover samples from Owl, assign 16 °C/0 °C groups, attach dissection covariates, write methylseq + RNA-seq sample sheets and the paired-sample set. |
| 26 | [`26-blood-RNAseq-align-DESeq2.Rmd`](26-blood-RNAseq-align-DESeq2.Rmd) | RNA-seq: trim → HISAT2 align to `GCF_031168955.1` → featureCounts → DESeq2 DEGs (16 vs 0 °C). |
| 27 | [`27-blood-WGBS-methylation-DMR.Rmd`](27-blood-WGBS-methylation-DMR.Rmd) | WGBS: nf-core/methylseq (Bismark) → CpG coverage → DMCs/DMRs (methylKit) → per-gene (gene-body + promoter) methylation. |
| 28 | [`28-blood-integration-feature-matrices.Rmd`](28-blood-integration-feature-matrices.Rmd) | Build matched expression × methylation matrices on shared genes and paired samples. |
| 29 | [`29-blood-integration-analysis.Rmd`](29-blood-integration-analysis.Rmd) | Integration: genome-wide meth↔expr relationship, per-gene correlation, **DEG × DMG overlap (candidate epigenetically-controlled genes)**, joint PCA / DIABLO. |
| 30 | [`30-blood-integration-figures.Rmd`](30-blood-integration-figures.Rmd) | Publication figures: genome-wide curve, joint PCA, correlation volcano, DEG×DMG quadrant, integration heatmap. |
| 31 | [`31-blood-integration-GO-enrichment.Rmd`](31-blood-integration-GO-enrichment.Rmd) | GO/functional enrichment of the DEG×DMG candidate set (and DEG/DMG sets), reusing the `03.2` gene→GO database; in-R topGO + DAVID export. |

### Notes

- This is a **separate track** from the earlier liver/gill/spleen RNA-seq (`05`–`13`) and
  methylation landscape (`18`, `20`) work; do not mix the blood count matrices with those.
- Heavy chunks default to `eval = FALSE` and assume UW Hyak (`klone`) paths — run them
  manually after staging reads and references off-repo (see the download chunk in `25`).
- The methylseq infrastructure (`16.2-genome-prep-klone`, `17.nextflow`, `17.config`) is
  reused for the blood WGBS run in `27`.
- Open question to confirm before trusting results: that `{NN}B` IDs map to the same fish
  across both assays, and that the 1–36/44–74 group assignment matches dissection records
  (notebook `25` cross-checks group vs recorded experiment temperature).
