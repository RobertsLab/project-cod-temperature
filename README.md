# project-cod-temperature

### Gene activity and genetic selection in Pacific cod reared under thermal stress

Juvenile Pacific cod will be reared in three temperatures under feeding and non-feeding conditions, then an integrated genomic approach will identify genes, gene variants, and epigenetic markers that respond to thermal stress and confer resilience. To complement the genomic approaches and further investigate temperature influences on energy resources, we will perform lipid analyses. This work will inform predictions of genetic selection and molecular response of Pacific cod in the Gulf of Alaska under climate change.

### Data

-   Experimental sampling data - [data/temp-experiment.csv](https://github.com/RobertsLab/project-cod-temperature/blob/main/data/temp-experiment.csv)

-   [WGS Sequencing Report](https://htmlpreview.github.io/?https://github.com/RobertsLab/project-cod-temperature/blob/main/output/Report_X202SC23041287-Z01-F001_20230611235113-4/X202SC23041287-Z01-F001_Report.html)

-   [WGS Raw data](https://owl.fish.washington.edu/nightingales/G_macrocephalus/H202SC23041287/01.RawData/)

-   [RNA sequencing raw reads](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/)

### How to work in this repo

The three main working directories are:

```         
data
code
output
```

For any code document, the name should start with a 2 number prefix (eg 01-temp-size-analysis.Rmd). All output from that code should be in a sub-directory of output named the same as the code. For example the output of 01-temp-size-analysis.Rmd would be in output/01-temp-size-analysis.Rmd/

### Analysis of RNA sequencing data

The workflow for RNAseq data analysis uses the following steps.

1.  Perform quality control (QC) and trimming on the raw reads (05-cod-RNAseq-trimming)
2.  Align trimmed reads to reference transcriptome and generate an estimate of transcript abundance in the form of a gene-level counts matrix (06-cod-RNAseq-alignment)
3.  Identify differentially expressed genes (DEGs), and generate associated visualizations (heatmap, volcano plot) (07-cod-RNAseq-DESeq2)
4.  Annotate reference transcriptome to generate a database of transcript IDs and associated gene ontology (GO) terms (03-transcriptome-annotation)
5.  Identify enriched GO terms (08-cod-RNAseq-GO-annotation)

### See also

<https://github.com/laurahspencer/pcod-juveniles-2023>
