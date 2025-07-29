# project-cod-temperature

### Gene activity and genetic selection in Pacific cod reared under thermal stress

Juvenile Pacific cod will be reared in three temperatures under feeding and non-feeding conditions, then an integrated genomic approach will identify genes, gene variants, and epigenetic markers that respond to thermal stress and confer resilience. To complement the genomic approaches and further investigate temperature influences on energy resources, we will perform lipid analyses. This work will inform predictions of genetic selection and molecular response of Pacific cod in the Gulf of Alaska under climate change.

### Data

- Heatwave juvenile cod rearing experiment data - [metadata-MHW-sample_collection_data.csv](data/metadata-MHW-sample_collection_data.csv)

-   Temperature experiment sampling data - [data/temp-experiment.csv](https://github.com/RobertsLab/project-cod-temperature/blob/main/data/temp-experiment.csv)

-   [Liver lipid data](https://github.com/RobertsLab/project-cod-temperature/blob/main/data/Lipid%20class%20liver%20data_091324.xlsx) from Louise 

-   [WGS Sequencing Report](https://htmlpreview.github.io/?https://github.com/RobertsLab/project-cod-temperature/blob/main/output/Report_X202SC23041287-Z01-F001_20230611235113-4/X202SC23041287-Z01-F001_Report.html)

-   [WGS Raw data](https://owl.fish.washington.edu/nightingales/G_macrocephalus/H202SC23041287/01.RawData/)

-   [RNA sequencing raw reads](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/)

-   **2025 WGS** /volume1/web/nightingales/G_macrocephalus/ [30-1149633765](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-1149633765/) and [30-1149634506](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-1149634506/)

_30-1149633765: heatwave genetics study of juvenile cod spanning 2008 - 2023 + six experimental fish that needed resequencing for juvenile temperature study.
30-1149634506: pilot cod ecotype sequencing (shallow, deep)_



### How to work in this repo

The three main working directories are:

```         
data
code
output
```

For any code document, the name should start with a 2 number prefix (eg 01-temp-size-analysis.Rmd). All output from that code should be in a sub-directory of output named the same as the code. For example the output of 01-temp-size-analysis.Rmd would be in output/01-temp-size-analysis.Rmd/

### Analysis of growth, condition, and lipids

Fish (n=40 per treatment, N=160) were tagged prior to the experiment, and per-individual growth metrics were collected to calculate specific growth rates based on wet weights (SGRww) and standard length (SGRsl) during experimental treatments, and Fulton's condition index based on wet weight (Kwet), and hepatosomatic index (HSI) at treatment termination. Please see [0-Phenotypes notebook](https://htmlpreview.github.io/?https://github.com/RobertsLab/project-cod-temperature/blob/main/general-notebooks/0-Phenotypes.html) for details. 

### Analysis of RNA sequencing data

The workflow for RNAseq data analysis uses the following steps. Visualizations of results can be viewed in the rendered `.md` files (e.g., 07-cod-RNAseq-DESeq2.md`)

1.  Perform quality control (QC) and trimming on the raw reads (`05-cod-RNAseq-trimming`). QC was performed using FastQC/MultiQC, and reads were trimmed using Flexbar.
2.  Align trimmed reads to reference transcriptome\genome and generate an estimate of transcript\gene abundance in the form of a gene-level counts matrix (`06-cod-RNAseq-alignment`, `06.2-cod-RNAseq-alignment-genome`). Reads were pseudoaligned to a transcriptome using kallisto, and transcript abundances were summarized to gene-level counts using Trinity abundance_estimates_to_matrix.pl. We also aligned reads to a genome using Hisat2 and summarized to gene-level counts using featureCounts.
3.  Identify differentially expressed genes (DEGs), and generate associated visualizations (heatmap, volcano plot) (`07-cod-RNAseq-DESeq2`, `07.2.1-cod-RNAseq-DESeq2genome-exon`, `07.2.2-cod-RNAseq-DESeq2-genome-gene`). Differential expression analysis was performed with DESeq2.
4.  Annotate reference transcriptome/genome to generate a database of transcript/gene IDs and associated gene ontology (GO) terms (`03-transcriptome-annotation`, `03.2-genome-annotation`)
5.  Identify enriched GO terms (`08-cod-RNAseq-GO-annotation`, `08.2.1-cod-RNAseq-GO-annotation-genome-exon`). Enrichment analysis was performed using DAVID. 

### See also

<https://github.com/laurahspencer/pcod-juveniles-2023>
