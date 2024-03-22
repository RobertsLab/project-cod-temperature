07-cod-RNAseq-DESeq2
================
Kathleen Durkin
2024-03-19

- <a href="#001-install-and-load-packages"
  id="toc-001-install-and-load-packages">0.0.1 Install and load
  packages</a>
- <a href="#1-deg-analysis-with-deseq2"
  id="toc-1-deg-analysis-with-deseq2">1 DEG Analysis with DESeq2</a>
  - <a href="#11-load-data" id="toc-11-load-data">1.1 Load data</a>
    - <a href="#111-load-count-data" id="toc-111-load-count-data">1.1.1 Load
      count data</a>
    - <a href="#112-import-sample-metadata-sheet"
      id="toc-112-import-sample-metadata-sheet">1.1.2 Import sample metadata
      sheet</a>
  - <a href="#12-preliminary-pca-visualization-liver-tissue"
    id="toc-12-preliminary-pca-visualization-liver-tissue">1.2 Preliminary
    PCA visualization (liver tissue)</a>
    - <a href="#121-deseq-object" id="toc-121-deseq-object">1.2.1 DESeq
      object</a>
    - <a href="#122-pca-visualization" id="toc-122-pca-visualization">1.2.2
      PCA visualization</a>
  - <a href="#13-treatment-comparisons-liver-tissue"
    id="toc-13-treatment-comparisons-liver-tissue">1.3 Treatment comparisons
    (liver tissue)</a>

Differential gene expression analysis for [Pacific cod RNAseq
data](https://shedurkin.github.io/Roberts-LabNotebook/posts/projects/pacific_cod/2023_12_13_pacific_cod.html).

- Raw reads found
  [here](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/)
- Reads aligned to transcriptome downloaded from
  [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_031168955.1/),
  stored
  [here](https://owl.fish.washington.edu/halfshell/genomic-databank/GCF_031168955.1_ASM3116895v1_rna.fna)
  as a part of lab [genomic
  resources](https://robertslab.github.io/resources/Genomic-Resources/#gadus-macrocephalus-pacific-cod).

### 0.0.1 Install and load packages

``` r
## clear
rm(list=ls())

## Install Rtools directly from (https://cran.r-project.org/bin/windows/Rtools/), then install these on first run:
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("vsn")
# BiocManager::install("tidybulk")
# BiocManager::install("goseq")
# BiocManager::install("affycoretools")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("pcaExplorer")
# BiocManager::install("apeglm")
# BiocManager::install("PCAtools")


# List of packages we want to install (run every time)
load.lib<-c("DESeq2","edgeR","goseq","dplyr","GenomicFeatures","data.table","calibrate","affycoretools","data.table","vsn","tidybulk","ggplot2","cowplot","pheatmap","gplots","RColorBrewer","EnhancedVolcano","pcaExplorer","readxl","apeglm","ashr","tibble","plotly","sqldf","PCAtools","ggpubr","beepr","genefilter","ComplexHeatmap","circlize","scales", "tidyverse")

# Select only the packages that aren't currently installed (run every time)
install.lib <- load.lib[!load.lib %in% installed.packages()]

# And finally we install the missing packages, including their dependency.
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
# After the installation process completes, we load all packages.
sapply(load.lib,require,character=TRUE)
```

# 1 DEG Analysis with DESeq2

## 1.1 Load data

### 1.1.1 Load count data

Load in the count matrix we generated after kallisto pseudoalignment
using the Trinity abundance_estimates_to_matrix.pl script. We also need
to slightly reformat the count matrix to make all of the estimated
counts integers, as required for DESeq2.

``` r
# Read in counts data. This is a gene-level counts matrix generated from kallisto transcript abundances using Trinity
cod_counts_data_OG <- read_delim("../output/06-cod-RNAseq-alignment/kallisto/kallisto.isoform.counts.matrix") 
head(cod_counts_data_OG)
```

``` r
# We need to modify this data frame so that the row names are actually row names, instead of comprising the first column
# Create object for row names from the first column
cod_row_names <- cod_counts_data_OG$...1

# Remove first column now that it's values are stored
# Note we're creating a new object for the data instead of modifying the original imported data
cod_counts_data <- cod_counts_data_OG[,-1]

# Assign the row names to the data using our object of stored names
rownames(cod_counts_data) <- cod_row_names

# Additional formatting
# Round all estimated counts to integers
cod_counts_data <- round(cod_counts_data, digits = 0)

# Remove the "kallisto_quant_" portion of the column names, to leave just the sample names
colnames(cod_counts_data) <- sub("kallisto_quant_", "sample_", colnames(cod_counts_data))

# Reorder the coumns into alphabetical order (to make it easier to create an associated metadata spreadsheet)
cod_counts_data <- cod_counts_data[, order(colnames(cod_counts_data))]

head(cod_counts_data)
```

### 1.1.2 Import sample metadata sheet

``` r
# Read in the csv as a data frame
cod_sample_info_OG <- read.csv("~/project-cod-temperature/data/DESeq2_Sample_Information.csv")
head(cod_sample_info_OG)

# Again, we need to reformat so that the data in the first column becomes the row names
sample_row_names <- cod_sample_info_OG$sample_name
cod_sample_info <- cod_sample_info_OG[,-1]
rownames(cod_sample_info) <- sample_row_names
head(cod_sample_info)
```

``` r
# Remove sample 92 (for now)
cod_sample_info <- cod_sample_info[rownames(cod_sample_info) != "sample_92", ]
head(cod_sample_info)

# Check that the column names of our count data are in the same order as the row names of our sample info sheet
ncol(cod_counts_data)
nrow(cod_sample_info)

colnames(cod_counts_data)==rownames(cod_sample_info)
```

## 1.2 Preliminary PCA visualization (liver tissue)

### 1.2.1 DESeq object

``` r
# Filter data
infosub_L <- cod_sample_info %>% filter(tissue_type == "Liver")
countsub_L <- subset(cod_counts_data, select=row.names(infosub_L))

# Calculate DESeq object
dds_L <- DESeqDataSetFromMatrix(countData = countsub_L,
                              colData = infosub_L,
                              design = ~ temp_treatment)

dds_L <- DESeq(dds_L)
resultsNames(dds_L) # lists the coefficients
```

### 1.2.2 PCA visualization

``` r
pca_L <- plotPCA(vst(dds_L), intgroup = c("temp_treatment"), returnData=TRUE)

percentVar_L <- round(100*attr(pca_L, "percentVar")) #plot PCA of samples with all data

# Assign specific colors to each temperature treatment level
temp_colors <- c(
  "0" = "darkblue",
  "5" = "royalblue1",
  "9" = "green",
  "16" = "orangered") 

ggplot(pca_L, aes(PC1, PC2, color=temp_treatment)) + 
  geom_point(size=4, alpha = 5/10) +
  xlab(paste0("PC1: ",percentVar_L[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_L[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=temp_colors)+
  stat_ellipse()

ggsave("../output/07-cod-RNAseq-DESeq2/Liver_allTemp_PCA.png")
```

## 1.3 Treatment comparisons (liver tissue)

The 9\*C temperature treatment is effectively our “control,” as it
represents the ambient temperature which wild juvenile Pacific cod would
experience.

``` r
# liver tissue, temperatures 9 vs. 16 

# Filter data
infosub_L_9_16 <- cod_sample_info %>% filter(tissue_type == "Liver" & (temp_treatment == "9" | temp_treatment == "16"))
countsub_L_9_16 <- subset(cod_counts_data, select=row.names(infosub_L_9_16))

# Calculate DESeq object
dds_L_9_16 <- DESeqDataSetFromMatrix(countData = countsub_L_9_16,
                              colData = infosub_L_9_16,
                              design = ~ temp_treatment)

# Convert temp_treatment to a factor (if it's not already)
dds_L_9_16$temp_treatment <- as.factor(dds_L_9_16$temp_treatment)

dds_L_9_16 <- DESeq(dds_L_9_16)
resultsNames(dds_L_9_16) # lists the coefficients
```

``` r
# Filtering: keep genes that have at least 10 counts across 1/3 of the samples - https://support.bioconductor.org/p/110307/
keep <- rowSums(DESeq2::counts(dds_L_9_16) >= 10) >= ncol(countsub_L_9_16)/3
dds_L_9_16<- dds_L_9_16[keep,]

# Generate Contrasts
contrast_list_9_16    <- c("temp_treatment", "16", "9") # order is important: factor, treatment group, control
res_table_L_9_16        <- results(dds_L_9_16, contrast=contrast_list_9_16, alpha = 0.05)
res_table_L_9_16_norm <- lfcShrink(dds_L_9_16,
                              coef=2, 
                              type="normal") # lfcThreshold = 0.585)  # a lfc threshold of 1 = 2-fold change, 0.585 = 1.5-fold change
res_table_L_9_16_apeglm <- lfcShrink(dds_L_9_16,
                              coef=2, 
                              type="apeglm") # lfcThreshold = 0.585)  # a lfc threshold of 1 = 2-fold change, 0.585 = 1.5-fold change
res_table_L_9_16_ashr   <- lfcShrink(dds_L_9_16,
                              coef=2, 
                              type="ashr")
```
