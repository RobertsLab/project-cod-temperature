06-cod-RNAseq-alignment
================
Kathleen Durkin
2024-03-19

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-align-to-reference-transcriptome-kallisto-pseudoalignment"
  id="toc-2-align-to-reference-transcriptome-kallisto-pseudoalignment">2
  Align to reference transcriptome (Kallisto pseudoalignment)</a>
  - <a href="#21-retrieving-the-reference-transcriptome"
    id="toc-21-retrieving-the-reference-transcriptome">2.1 Retrieving the
    reference transcriptome</a>
  - <a href="#22-verify-transcriptome-fasta-md5-checksum"
    id="toc-22-verify-transcriptome-fasta-md5-checksum">2.2 Verify
    transcriptome FastA MD5 checksum</a>
  - <a href="#23-building-index" id="toc-23-building-index">2.3 Building
    Index</a>
  - <a href="#24-sample-quantification"
    id="toc-24-sample-quantification">2.4 Sample Quantification</a>
  - <a href="#25-multiqc-on-kallisto-output-logs"
    id="toc-25-multiqc-on-kallisto-output-logs">2.5 MultiQC on Kallisto
    output logs</a>
  - <a href="#26-trinity-matrix-with-kallisto-output"
    id="toc-26-trinity-matrix-with-kallisto-output">2.6 Trinity Matrix with
    Kallisto Output</a>

``` r
library(tidyverse)
library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(plotly)
```

Code for aligning RNAseq data to reference transcriptome/genome, to be
used on [Pacific cod RNAseq
data](https://shedurkin.github.io/Roberts-LabNotebook/posts/projects/pacific_cod/2023_12_13_pacific_cod.html).

- Raw reads found
  [here](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/)
- Transcriptome downloaded from
  [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_031168955.1/),
  stored
  [here](https://owl.fish.washington.edu/halfshell/genomic-databank/GCF_031168955.1_ASM3116895v1_rna.fna)
  as a part of lab [genomic
  resources](https://robertslab.github.io/resources/Genomic-Resources/#gadus-macrocephalus-pacific-cod)

Note: Kallisto pseudoalignment doesn’t necessarily require input reads
to be trimmed, provided they are of sufficient quality.

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export cod_dir=/home/shared/8TB_HDD_02/shedurkin/project-cod-temperature'
echo 'export output_dir_top=${cod_dir}/output/06-cod-RNAseq-alignment'
echo 'export raw_fastqc_dir=${cod_dir}/output/05-cod-RNAseq-trimming/raw-fastqc'
echo 'export raw_reads_dir=${cod_dir}/data/05-cod-RNAseq-trimming/raw-reads'
echo 'export trimmed_fastqc_dir=${cod_dir}/output/05-cod-RNAseq-trimming/trimmed-fastqc'
echo 'export trimmed_reads_dir=${cod_dir}/output/05-cod-RNAseq-trimming/trimmed-reads'
echo 'export kallisto_output_dir=${output_dir_top}/kallisto'
echo ""


echo "# Input/Output files"
echo 'export transcriptome_fasta_dir=${cod_dir}/data'
echo 'export transcriptome_fasta_name="GCF_031168955.1_ASM3116895v1_rna"'
echo 'export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"'
echo 'export kallisto_index_name="G_macrocephalus_kallisto_index.idx"'


echo "# External data URLs and checksums"
echo 'export transcriptome_fasta_url="https://owl.fish.washington.edu/halfshell/genomic-databank/GCF_031168955.1_ASM3116895v1_rna.fna"'
echo 'export transcriptome_checksum="2a6c7c98982727e688f033a9b236725b"'
echo ""


echo "# Paths to programs"
echo 'export kallisto=/home/shared/kallisto/kallisto'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo 'export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl'
echo ""


echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""


echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[kallisto]="${kallisto}" \'
echo '[multiqc]="${multiqc}" \'
echo '[trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \'
echo ")"
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export cod_dir=/home/shared/8TB_HDD_02/shedurkin/project-cod-temperature
    export output_dir_top=${cod_dir}/output/06-cod-RNAseq-alignment
    export raw_fastqc_dir=${cod_dir}/output/05-cod-RNAseq-trimming/raw-fastqc
    export raw_reads_dir=${cod_dir}/data/05-cod-RNAseq-trimming/raw-reads
    export trimmed_fastqc_dir=${cod_dir}/output/05-cod-RNAseq-trimming/trimmed-fastqc
    export trimmed_reads_dir=${cod_dir}/output/05-cod-RNAseq-trimming/trimmed-reads
    export kallisto_output_dir=${output_dir_top}/kallisto

    # Input/Output files
    export transcriptome_fasta_dir=${cod_dir}/data
    export transcriptome_fasta_name="GCF_031168955.1_ASM3116895v1_rna"
    export transcriptome_fasta="${transcriptome_fasta_dir}/${transcriptome_fasta_name}"
    export kallisto_index_name="G_macrocephalus_kallisto_index.idx"
    # External data URLs and checksums
    export transcriptome_fasta_url="https://owl.fish.washington.edu/halfshell/genomic-databank/GCF_031168955.1_ASM3116895v1_rna.fna"
    export transcriptome_checksum="2a6c7c98982727e688f033a9b236725b"

    # Paths to programs
    export kallisto=/home/shared/kallisto/kallisto
    export multiqc=/home/sam/programs/mambaforge/bin/multiqc
    export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl

    # Set number of CPUs to use
    export threads=20

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [kallisto]="${kallisto}" \
    [multiqc]="${multiqc}" \
    [trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \
    )

There didn’t seem to be a significant difference in sequence quality
following trimming, so for now I’m proceeding with the raw reads (though
I’ll likely eventually rerun the kallisto with trimmed)

# 2 Align to reference transcriptome (Kallisto pseudoalignment)

## 2.1 Retrieving the reference transcriptome

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${transcriptome_fasta_dir} \
--recursive \
--no-check-certificate \
--continue \
--no-host-directories \
--no-directories \
--no-parent \
--quiet \
--execute robots=off \
--accept "${transcriptome_fasta_name}.fna" ${transcriptome_fasta_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${transcriptome_fasta_dir}"
```

    total 1.9G
    drwxr-xr-x 3 shedurkin labmembers 4.0K Mar  4 11:05 05-cod-RNAseq-trimming
    -rw-r--r-- 1 shedurkin labmembers  13K Dec 27 15:45 Cod_RNAseq_NGS_Template_File.xlsx
    -rw-r--r-- 1 shedurkin labmembers 2.1K Mar 20 20:55 DESeq2_Sample_Information.csv
    -rw-r--r-- 1 shedurkin labmembers  38M Oct 25  2023 Gadus_macrocephalus.coding.gene.V1.cds
    -rw-r--r-- 1 shedurkin labmembers 537M Oct 16  2023 GCF_031168955.1_ASM3116895v1_genomic.fna
    -rw-r--r-- 1 shedurkin labmembers 351M Oct 16  2023 GCF_031168955.1_ASM3116895v1.gff
    -rw-r--r-- 1 shedurkin labmembers 169M Oct 16  2023 GCF_031168955.1_ASM3116895v1_rna.fna
    -rw-r--r-- 1 shedurkin labmembers 404M Apr 23 14:29 genomic.gtf
    -rw-r--r-- 1 shedurkin labmembers  47K Oct 25  2023 Pcod Temp Growth experiment 2022-23 DATA.xlsx
    -rw-r--r-- 1 shedurkin labmembers 231K Mar  4 17:41 Sample.QC.report.of_30-943133806_240118025106.pdf
    -rw-r--r-- 1 shedurkin labmembers  12K Mar  4 17:41 Sample.QC.report.of_30-943133806_240118025106.xlsx
    -rw-r--r-- 1 shedurkin labmembers  12K Oct 25  2023 temp-experiment.csv
    -rw-r--r-- 1 shedurkin labmembers 271M Oct 25  2023 uniprot_sprot_r2023_04.fasta
    -rw-r--r-- 1 shedurkin labmembers  88M Apr 17 11:54 uniprot_sprot_r2023_04.fasta.gz

## 2.2 Verify transcriptome FastA MD5 checksum

``` bash
# Load bash variables into memory
source .bashvars

cd "${transcriptome_fasta_dir}"

# Checksums file contains other files, so this just looks for the sRNAseq files.
md5sum --check <<< "${transcriptome_checksum}  ${transcriptome_fasta_name}.fna"
```

    GCF_031168955.1_ASM3116895v1_rna.fna: OK

## 2.3 Building Index

``` bash
# Load bash variables into memory
source .bashvars

cd "${kallisto_output_dir}"

${programs_array[kallisto]} index \
--threads=${threads} \
--index="${kallisto_index_name}" \
"${transcriptome_fasta}.fna"
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh ${kallisto_output_dir}
```

    total 1.5G
    -rw-r--r-- 1 shedurkin labmembers 1.5G Mar 18 16:08 G_macrocephalus_kallisto_index.idx
    -rw-r--r-- 1 shedurkin labmembers  20M Mar 19 16:45 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Mar 19 16:45 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers  24M Mar 19 16:45 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Mar 19 16:45 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:46 kallisto_quant_1
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:11 kallisto_quant_10
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:03 kallisto_quant_100
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:04 kallisto_quant_100.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:05 kallisto_quant_107
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:05 kallisto_quant_107.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:07 kallisto_quant_108
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:08 kallisto_quant_108.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:09 kallisto_quant_109
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:10 kallisto_quant_109.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:12 kallisto_quant_10.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:22 kallisto_quant_11
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:13 kallisto_quant_110
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:14 kallisto_quant_110.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:15 kallisto_quant_117
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:16 kallisto_quant_117.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:18 kallisto_quant_118
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:18 kallisto_quant_118.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:20 kallisto_quant_119
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:20 kallisto_quant_119.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:22 kallisto_quant_11.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:33 kallisto_quant_12
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:24 kallisto_quant_120
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:24 kallisto_quant_120.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:26 kallisto_quant_121
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:27 kallisto_quant_121.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:28 kallisto_quant_127
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:29 kallisto_quant_127.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:30 kallisto_quant_128
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:30 kallisto_quant_128.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:31 kallisto_quant_129
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:31 kallisto_quant_129.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:33 kallisto_quant_12.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:39 kallisto_quant_13
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:34 kallisto_quant_131
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:34 kallisto_quant_131.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:35 kallisto_quant_137
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:35 kallisto_quant_137.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:36 kallisto_quant_138
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:36 kallisto_quant_138.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:37 kallisto_quant_139
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:37 kallisto_quant_139.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:39 kallisto_quant_13.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:40 kallisto_quant_140
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:40 kallisto_quant_140.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:41 kallisto_quant_147
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:41 kallisto_quant_147.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:43 kallisto_quant_148
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:44 kallisto_quant_148.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:51 kallisto_quant_149
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:52 kallisto_quant_149.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:52 kallisto_quant_150
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 22:53 kallisto_quant_150.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:37 kallisto_quant_18
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:37 kallisto_quant_18.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:41 kallisto_quant_19
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:39 kallisto_quant_19-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:40 kallisto_quant_19-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:41 kallisto_quant_19.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:43 kallisto_quant_19-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:44 kallisto_quant_19-S.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:46 kallisto_quant_1.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:59 kallisto_quant_2
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:50 kallisto_quant_20
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:48 kallisto_quant_20-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:49 kallisto_quant_20-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:50 kallisto_quant_20.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:52 kallisto_quant_20-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:52 kallisto_quant_20-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:54 kallisto_quant_21
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:54 kallisto_quant_21.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:56 kallisto_quant_28
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:56 kallisto_quant_28.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:57 kallisto_quant_29
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:58 kallisto_quant_29.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:00 kallisto_quant_2.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:53 kallisto_quant_3
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 18:01 kallisto_quant_30
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:02 kallisto_quant_30.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 18:03 kallisto_quant_31
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:04 kallisto_quant_31.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:47 kallisto_quant_37
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:48 kallisto_quant_37.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:49 kallisto_quant_38
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:49 kallisto_quant_38.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:51 kallisto_quant_39
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:51 kallisto_quant_39.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:53 kallisto_quant_3.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:04 kallisto_quant_4
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:54 kallisto_quant_40
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:55 kallisto_quant_40.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:56 kallisto_quant_41
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:56 kallisto_quant_41.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:58 kallisto_quant_47
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:58 kallisto_quant_47.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:00 kallisto_quant_48
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:00 kallisto_quant_48.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:02 kallisto_quant_49
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:02 kallisto_quant_49.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:04 kallisto_quant_4.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:26 kallisto_quant_5
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:06 kallisto_quant_50
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:06 kallisto_quant_50.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:14 kallisto_quant_57
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:10 kallisto_quant_57-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:11 kallisto_quant_57-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:14 kallisto_quant_57.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:15 kallisto_quant_57-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:15 kallisto_quant_57-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:20 kallisto_quant_58
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:18 kallisto_quant_58-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:18 kallisto_quant_58-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:20 kallisto_quant_58.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:22 kallisto_quant_58-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:22 kallisto_quant_58-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:24 kallisto_quant_59
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:24 kallisto_quant_59.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:26 kallisto_quant_5.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:28 kallisto_quant_60
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:28 kallisto_quant_60.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:30 kallisto_quant_67
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:30 kallisto_quant_67.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:32 kallisto_quant_68
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:32 kallisto_quant_68.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:34 kallisto_quant_69
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:34 kallisto_quant_69.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:35 kallisto_quant_70
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:36 kallisto_quant_70.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:37 kallisto_quant_78
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:37 kallisto_quant_78.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:39 kallisto_quant_79
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:39 kallisto_quant_79.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:40 kallisto_quant_80
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:41 kallisto_quant_80.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:42 kallisto_quant_83
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:43 kallisto_quant_83.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:44 kallisto_quant_88
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:44 kallisto_quant_88.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:46 kallisto_quant_90
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:46 kallisto_quant_90.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:48 kallisto_quant_91
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:48 kallisto_quant_91.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:50 kallisto_quant_97
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:50 kallisto_quant_97.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:52 kallisto_quant_98
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:52 kallisto_quant_98.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:53 kallisto_quant_99
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:54 kallisto_quant_99.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:55 kallisto_quant_RESUB-116
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:55 kallisto_quant_RESUB-116.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:57 kallisto_quant_RESUB-156
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:57 kallisto_quant_RESUB-156.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:58 kallisto_quant_RESUB-36
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:59 kallisto_quant_RESUB-36.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:54 kallisto_quant_RESUB-76
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 22:54 kallisto_quant_RESUB-76.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:13 kallisto_quant_RESUB-94
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 22:14 kallisto_quant_RESUB-94.log

## 2.4 Sample Quantification

Kallisto can run quantification on either single- or paired-end reads.
The default option is paired-end, which requires the input of an even
number of paired fastq files (e.g., pairA_R1.fastq, pairA_R2.fastq). To
use single-end mode, include the –single flag, as well as -l
(–fragment-length=DOUBLE, estimated avg. fragment length) and -s
(–sd=DOUBLE, estimates stand. dev. of fragment length), and a number of
fastq files. Again, gzipped files are acceptable.

Kallisto quant is rather finicky about how you input sets of paired
reads, and you can only input a single pair at a time. To circumvent,
I’ll create a quantification function and apply it iteratively to each
pair of reads using a loop.

``` bash
# Load bash variables into memory
source .bashvars

# Function to run kallisto quant. Takes two (paired) reads as input, outputs to sample-associated directory
run_kallisto_quant() {
    source .bashvars  # Source .bashvars inside the function to make its variables accessible
    local R1_fastq=${1}
    local R2_fastq=${2}
    
    cd ${kallisto_output_dir}
    sample_num=$(basename "${R1_fastq}" "_R1_001.fastq.gz")
    mkdir kallisto_quant_${sample_num}

    ${programs_array[kallisto]} quant \
        --threads=${threads} \
        --index="${kallisto_output_dir}/${kallisto_index_name}" \
        --output-dir="${kallisto_output_dir}/kallisto_quant_${sample_num}" \
        --bootstrap-samples=100 \
        ${raw_reads_dir}/${R1_fastq} ${raw_reads_dir}/${R2_fastq} \
        &> "${kallisto_output_dir}/kallisto_quant_${sample_num}.log"
}



# Iteratively apply run_kallisto_quant on each pair of input reads
for file_r1 in "${raw_reads_dir}"/*_R1_001.fastq.gz; do
    # Extract the sample name from the file name
    sample_name=$(basename "${file_r1}" "_R1_001.fastq.gz")

    # Form the file names (function takes input file names, not paths)
    file_r1_name="${sample_name}_R1_001.fastq.gz"
    file_r2_name="${sample_name}_R2_001.fastq.gz"

    # Check that the sample hasn't already been quantified
    if [ ! -d "${kallisto_output_dir}/kallisto_quant_${sample_name}" ]; then
    
        # Check if the corresponding R2 file exists
        if [ -e "${raw_reads_dir}/${file_r2}" ]; then
            # Run kallisto quant on the file pair
            run_kallisto_quant "${file_r1_name}" "${file_r2_name}" 

            echo "Processed sample: ${sample_name}"
        fi
    else
        echo "Sample already processed: ${sample_name}"
    fi
done
```

Check that we have the appropriate number of output folders. We should
have one log file for each pair of reads

``` bash
# Load bash variables into memory
source .bashvars

# Count number of raw read files
cd ${raw_reads_dir}
echo "Number of raw reads:"
ls -1 | wc -l

# Count number of kallisto output 
cd ${kallisto_output_dir}
echo "Number of output log files"
find . -type f -name "*.log" | wc -l
```

    Number of raw reads:
    158
    Number of output log files
    79

## 2.5 MultiQC on Kallisto output logs

``` bash
# Load bash variables into memory
source .bashvars

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${kallisto_output_dir}/*.log -o ${output_dir_top}

echo ""
echo "MultiQC on raw FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${output_dir_top}/*.zip
echo "FastQC zip files removed."
echo ""

# View directory contents
ls -lh ${output_dir_top}
```

    Beginning MultiQC on raw FastQC...


      /// MultiQC 🔍 | v1.14

    |           multiqc | MultiQC Version v1.21 now available!
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_100.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_107.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_108.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_109.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_10.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_110.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_117.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_118.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_119.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_11.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_120.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_121.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_127.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_128.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_129.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_12.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_131.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_137.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_138.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_139.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_13.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_140.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_147.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_148.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_149.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_150.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_18.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19-G.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19-S.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_1.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20-G.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20-S.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_21.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_28.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_29.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_2.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_30.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_31.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_37.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_38.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_39.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_3.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_40.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_41.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_47.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_48.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_49.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_4.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_50.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57-G.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57-S.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58-G.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58-S.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_59.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_5.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_60.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_67.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_68.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_69.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_70.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_78.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_79.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_80.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_83.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_88.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_90.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_91.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_97.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_98.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_99.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-116.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-156.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-36.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-76.log
    |           multiqc | Search path : /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-94.log
    |         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 79/79  
    |          kallisto | Found 79 reports
    |           multiqc | Compressing plot data
    |           multiqc | Previous MultiQC output found! Adjusting filenames..
    |           multiqc | Use -f or --force to overwrite existing reports instead
    |           multiqc | Report      : ../output/06-cod-RNAseq-alignment/multiqc_report_1.html
    |           multiqc | Data        : ../output/06-cod-RNAseq-alignment/multiqc_data_1
    |           multiqc | MultiQC complete

    MultiQC on raw FastQs complete.

    Removing FastQC zip files.

    rm: cannot remove '/home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/*.zip': No such file or directory
    FastQC zip files removed.

    total 2.3M
    drwxr-xr-x 81 shedurkin labmembers  12K Apr 29 22:53 kallisto
    drwxr-xr-x  2 shedurkin labmembers 4.0K Apr 29 22:55 multiqc_data
    drwxr-xr-x  2 shedurkin labmembers 4.0K Apr 29 23:20 multiqc_data_1
    -rw-r--r--  1 shedurkin labmembers 1.2M Apr 29 23:20 multiqc_report_1.html
    -rw-r--r--  1 shedurkin labmembers 1.2M Apr 29 22:55 multiqc_report.html

I also want to include the treatment/tank info when plotting alignment
rates across samples

``` r
# Load multiqc stats
kallisto_multiqc <- read.csv("../output/06-cod-RNAseq-alignment/multiqc_data/multiqc_kallisto.txt", sep = '\t')
# Adjust sample name formatting (to prep for join)
kallisto_multiqc$Sample <- gsub("_R1_001", "", kallisto_multiqc$Sample) 
kallisto_multiqc$Sample <- paste("sample_", kallisto_multiqc$Sample, sep = "")
# Load experimental data
cod_sample_info_OG <- read.csv("../data/DESeq2_Sample_Information.csv")

kallisto_multiqc_plustreatment <- left_join(cod_sample_info_OG, kallisto_multiqc, by = c("sample_name" = "Sample")) %>% 
  na.omit()
kallisto_multiqc_plustreatment <- kallisto_multiqc_plustreatment[order(kallisto_multiqc_plustreatment$sample_number),]

a <- ggplot(kallisto_multiqc_plustreatment,
       aes(x=reorder(sample_name, sample_number), y=percent_aligned, fill=as.factor(temp_treatment))) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=7))

ggplotly(a)
```

<div class="plotly html-widget html-fill-item" id="htmlwidget-d77b47e2a16e7707e9f9" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-d77b47e2a16e7707e9f9">{"x":{"data":[{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.90000000000000213,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],"y":[78.137666740364068,72.289551451902852,75.65512701806442,70.292892735598585,63.887216917731735,76.052982059048475,75.920547971626178,72.836534408040066,74.506506231057386,78.373689421732564,75.29179932210657,39.498585386987685,65.78697556179219,74.564147342965597,68.48448055551313,79.012048738854773,73.069596027789444,75.229249487301502,72.720421898167999,74.399464671647593,73.692177330834483,75.472746916236488],"text":["reorder(sample_name, sample_number): sample_37<br />percent_aligned: 78.13767<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_38<br />percent_aligned: 72.28955<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_39<br />percent_aligned: 75.65513<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_40<br />percent_aligned: 70.29289<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_41<br />percent_aligned: 63.88722<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_47<br />percent_aligned: 76.05298<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_48<br />percent_aligned: 75.92055<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_49<br />percent_aligned: 72.83653<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_50<br />percent_aligned: 74.50651<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_57<br />percent_aligned: 78.37369<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_57-G<br />percent_aligned: 75.29180<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_57-S<br />percent_aligned: 39.49859<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_58<br />percent_aligned: 65.78698<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_58-G<br />percent_aligned: 74.56415<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_58-S<br />percent_aligned: 68.48448<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_59<br />percent_aligned: 79.01205<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_60<br />percent_aligned: 73.06960<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_67<br />percent_aligned: 75.22925<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_68<br />percent_aligned: 72.72042<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_69<br />percent_aligned: 74.39946<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_70<br />percent_aligned: 73.69218<br />as.factor(temp_treatment): 0","reorder(sample_name, sample_number): sample_RESUB-76<br />percent_aligned: 75.47275<br />as.factor(temp_treatment): 0"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"0","legendgroup":"0","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79],"y":[77.021005375035699,80.757990130084593,83.778896934841157,78.350538434832899,72.406726773047595,77.201901733727951,74.27127985486031,37.513100408341273,71.607759191415042,75.850130921805331,77.935136758050177,78.109721850451876,75.895593344802421,74.739847670978776,71.385380896714764,58.313748092544493,69.230076755631458,79.324083185693055],"text":["reorder(sample_name, sample_number): sample_117<br />percent_aligned: 77.02101<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_118<br />percent_aligned: 80.75799<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_119<br />percent_aligned: 83.77890<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_120<br />percent_aligned: 78.35054<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_121<br />percent_aligned: 72.40673<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_127<br />percent_aligned: 77.20190<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_128<br />percent_aligned: 74.27128<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_129<br />percent_aligned: 37.51310<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_131<br />percent_aligned: 71.60776<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_137<br />percent_aligned: 75.85013<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_138<br />percent_aligned: 77.93514<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_139<br />percent_aligned: 78.10972<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_140<br />percent_aligned: 75.89559<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_147<br />percent_aligned: 74.73985<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_148<br />percent_aligned: 71.38538<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_149<br />percent_aligned: 58.31375<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_150<br />percent_aligned: 69.23008<br />as.factor(temp_treatment): 5","reorder(sample_name, sample_number): sample_RESUB-156<br />percent_aligned: 79.32408<br />as.factor(temp_treatment): 5"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(124,174,0,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"5","legendgroup":"5","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61],"y":[76.078769568721725,79.148045486106852,74.173493451808454,79.228327961793738,63.335167636042634,70.964840711590696,75.422840105451399,78.309635680699046,78.01720630789174,69.647810861164388,71.94534440733608,79.27637489177765,79.166188449495436,68.727363448190431,78.212098206040622,81.473563533017696,77.510652847664048],"text":["reorder(sample_name, sample_number): sample_78<br />percent_aligned: 76.07877<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_79<br />percent_aligned: 79.14805<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_80<br />percent_aligned: 74.17349<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_83<br />percent_aligned: 79.22833<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_88<br />percent_aligned: 63.33517<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_90<br />percent_aligned: 70.96484<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_91<br />percent_aligned: 75.42284<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_RESUB-94<br />percent_aligned: 78.30964<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_97<br />percent_aligned: 78.01721<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_98<br />percent_aligned: 69.64781<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_99<br />percent_aligned: 71.94534<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_100<br />percent_aligned: 79.27637<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_107<br />percent_aligned: 79.16619<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_108<br />percent_aligned: 68.72736<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_109<br />percent_aligned: 78.21210<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_110<br />percent_aligned: 81.47356<br />as.factor(temp_treatment): 9","reorder(sample_name, sample_number): sample_RESUB-116<br />percent_aligned: 77.51065<br />as.factor(temp_treatment): 9"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"9","legendgroup":"9","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999947,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],"y":[76.984509009072383,81.217597602966492,72.710852216810622,75.270654754232197,75.525978232364537,73.086801075491749,78.684582297934085,75.320312039981246,70.729288335932111,72.693136519714258,76.974505812254847,70.67992703598685,74.274929233177062,71.359779859013699,70.735886228605153,73.7971666088321,73.544591341789584,76.298847432170888,69.238566322686296,74.978069400112346,80.026266199246891,82.193330042352528],"text":["reorder(sample_name, sample_number): sample_1<br />percent_aligned: 76.98451<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_2<br />percent_aligned: 81.21760<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_3<br />percent_aligned: 72.71085<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_4<br />percent_aligned: 75.27065<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_5<br />percent_aligned: 75.52598<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_10<br />percent_aligned: 73.08680<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_11<br />percent_aligned: 78.68458<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_12<br />percent_aligned: 75.32031<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_13<br />percent_aligned: 70.72929<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_18<br />percent_aligned: 72.69314<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_19<br />percent_aligned: 76.97451<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_19-G<br />percent_aligned: 70.67993<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_19-S<br />percent_aligned: 74.27493<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_20<br />percent_aligned: 71.35978<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_20-G<br />percent_aligned: 70.73589<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_20-S<br />percent_aligned: 73.79717<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_21<br />percent_aligned: 73.54459<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_28<br />percent_aligned: 76.29885<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_29<br />percent_aligned: 69.23857<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_30<br />percent_aligned: 74.97807<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_31<br />percent_aligned: 80.02627<br />as.factor(temp_treatment): 16","reorder(sample_name, sample_number): sample_RESUB-36<br />percent_aligned: 82.19333<br />as.factor(temp_treatment): 16"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(199,124,255,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"16","legendgroup":"16","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.228310502283104,"r":7.3059360730593621,"b":40.05685961550482,"l":37.260273972602747},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,79.599999999999994],"tickmode":"array","ticktext":["sample_1","sample_2","sample_3","sample_4","sample_5","sample_10","sample_11","sample_12","sample_13","sample_18","sample_19","sample_19-G","sample_19-S","sample_20","sample_20-G","sample_20-S","sample_21","sample_28","sample_29","sample_30","sample_31","sample_RESUB-36","sample_37","sample_38","sample_39","sample_40","sample_41","sample_47","sample_48","sample_49","sample_50","sample_57","sample_57-G","sample_57-S","sample_58","sample_58-G","sample_58-S","sample_59","sample_60","sample_67","sample_68","sample_69","sample_70","sample_RESUB-76","sample_78","sample_79","sample_80","sample_83","sample_88","sample_90","sample_91","sample_RESUB-94","sample_97","sample_98","sample_99","sample_100","sample_107","sample_108","sample_109","sample_110","sample_RESUB-116","sample_117","sample_118","sample_119","sample_120","sample_121","sample_127","sample_128","sample_129","sample_131","sample_137","sample_138","sample_139","sample_140","sample_147","sample_148","sample_149","sample_150","sample_RESUB-156"],"tickvals":[1,2,3,4,5,6,7,8,9,10,10.999999999999998,12,13,14,15,16,17,18,19,20,21,22,23,24.000000000000004,25,26,27,28,29,30,31.000000000000004,32,33,34,35,36,37,38,39,40,41.000000000000007,42,43,44,45,46.000000000000007,47,48,49,50,50.999999999999993,52,53,54,55.000000000000007,56,57,58,59,59.999999999999993,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79],"categoryorder":"array","categoryarray":["sample_1","sample_2","sample_3","sample_4","sample_5","sample_10","sample_11","sample_12","sample_13","sample_18","sample_19","sample_19-G","sample_19-S","sample_20","sample_20-G","sample_20-S","sample_21","sample_28","sample_29","sample_30","sample_31","sample_RESUB-36","sample_37","sample_38","sample_39","sample_40","sample_41","sample_47","sample_48","sample_49","sample_50","sample_57","sample_57-G","sample_57-S","sample_58","sample_58-G","sample_58-S","sample_59","sample_60","sample_67","sample_68","sample_69","sample_70","sample_RESUB-76","sample_78","sample_79","sample_80","sample_83","sample_88","sample_90","sample_91","sample_RESUB-94","sample_97","sample_98","sample_99","sample_100","sample_107","sample_108","sample_109","sample_110","sample_RESUB-116","sample_117","sample_118","sample_119","sample_120","sample_121","sample_127","sample_128","sample_129","sample_131","sample_137","sample_138","sample_139","sample_140","sample_147","sample_148","sample_149","sample_150","sample_RESUB-156"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":9.2984640929846396},"tickangle":-60,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"reorder(sample_name, sample_number)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-4.188944846742058,87.967841781583218],"tickmode":"array","ticktext":["0","20","40","60","80"],"tickvals":[0,20,40,60,80],"categoryorder":"array","categoryarray":["0","20","40","60","80"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"percent_aligned","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"title":{"text":"as.factor(temp_treatment)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"a1a895490de":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"a1a895490de","visdat":{"a1a895490de":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

``` r
b <- ggplot(kallisto_multiqc_plustreatment,
       aes(x=reorder(sample_name, sample_number), y=percent_aligned, fill=as.factor(tank))) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=7))

ggplotly(b)
```

<div class="plotly html-widget html-fill-item" id="htmlwidget-10760e25d67332378df1" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-10760e25d67332378df1">{"x":{"data":[{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0,0],"x":[1,2,3,4,5],"y":[76.984509009072383,81.217597602966492,72.710852216810622,75.270654754232197,75.525978232364537],"text":["reorder(sample_name, sample_number): sample_1<br />percent_aligned: 76.98451<br />as.factor(tank): 1","reorder(sample_name, sample_number): sample_2<br />percent_aligned: 81.21760<br />as.factor(tank): 1","reorder(sample_name, sample_number): sample_3<br />percent_aligned: 72.71085<br />as.factor(tank): 1","reorder(sample_name, sample_number): sample_4<br />percent_aligned: 75.27065<br />as.factor(tank): 1","reorder(sample_name, sample_number): sample_5<br />percent_aligned: 75.52598<br />as.factor(tank): 1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"1","legendgroup":"1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.89999999999999947,0.89999999999999858],"base":[0,0,0,0],"x":[6,7,8,9],"y":[73.086801075491749,78.684582297934085,75.320312039981246,70.729288335932111],"text":["reorder(sample_name, sample_number): sample_10<br />percent_aligned: 73.08680<br />as.factor(tank): 2","reorder(sample_name, sample_number): sample_11<br />percent_aligned: 78.68458<br />as.factor(tank): 2","reorder(sample_name, sample_number): sample_12<br />percent_aligned: 75.32031<br />as.factor(tank): 2","reorder(sample_name, sample_number): sample_13<br />percent_aligned: 70.72929<br />as.factor(tank): 2"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(230,134,19,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"2","legendgroup":"2","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0,0,0,0],"x":[10,11,12,13,14,15,16,17],"y":[72.693136519714258,76.974505812254847,70.67992703598685,74.274929233177062,71.359779859013699,70.735886228605153,73.7971666088321,73.544591341789584],"text":["reorder(sample_name, sample_number): sample_18<br />percent_aligned: 72.69314<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_19<br />percent_aligned: 76.97451<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_19-G<br />percent_aligned: 70.67993<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_19-S<br />percent_aligned: 74.27493<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_20<br />percent_aligned: 71.35978<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_20-G<br />percent_aligned: 70.73589<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_20-S<br />percent_aligned: 73.79717<br />as.factor(tank): 3","reorder(sample_name, sample_number): sample_21<br />percent_aligned: 73.54459<br />as.factor(tank): 3"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(205,150,0,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"3","legendgroup":"3","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0],"x":[18,19,20,21,22],"y":[76.298847432170888,69.238566322686296,74.978069400112346,80.026266199246891,82.193330042352528],"text":["reorder(sample_name, sample_number): sample_28<br />percent_aligned: 76.29885<br />as.factor(tank): 4","reorder(sample_name, sample_number): sample_29<br />percent_aligned: 69.23857<br />as.factor(tank): 4","reorder(sample_name, sample_number): sample_30<br />percent_aligned: 74.97807<br />as.factor(tank): 4","reorder(sample_name, sample_number): sample_31<br />percent_aligned: 80.02627<br />as.factor(tank): 4","reorder(sample_name, sample_number): sample_RESUB-36<br />percent_aligned: 82.19333<br />as.factor(tank): 4"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(171,163,0,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"4","legendgroup":"4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0],"x":[62,63,64,65,66],"y":[77.021005375035699,80.757990130084593,83.778896934841157,78.350538434832899,72.406726773047595],"text":["reorder(sample_name, sample_number): sample_117<br />percent_aligned: 77.02101<br />as.factor(tank): 5","reorder(sample_name, sample_number): sample_118<br />percent_aligned: 80.75799<br />as.factor(tank): 5","reorder(sample_name, sample_number): sample_119<br />percent_aligned: 83.77890<br />as.factor(tank): 5","reorder(sample_name, sample_number): sample_120<br />percent_aligned: 78.35054<br />as.factor(tank): 5","reorder(sample_name, sample_number): sample_121<br />percent_aligned: 72.40673<br />as.factor(tank): 5"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(124,174,0,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"5","legendgroup":"5","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0],"x":[67,68,69,70],"y":[77.201901733727951,74.27127985486031,37.513100408341273,71.607759191415042],"text":["reorder(sample_name, sample_number): sample_127<br />percent_aligned: 77.20190<br />as.factor(tank): 6","reorder(sample_name, sample_number): sample_128<br />percent_aligned: 74.27128<br />as.factor(tank): 6","reorder(sample_name, sample_number): sample_129<br />percent_aligned: 37.51310<br />as.factor(tank): 6","reorder(sample_name, sample_number): sample_131<br />percent_aligned: 71.60776<br />as.factor(tank): 6"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(12,183,2,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"6","legendgroup":"6","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0],"x":[71,72,73,74],"y":[75.850130921805331,77.935136758050177,78.109721850451876,75.895593344802421],"text":["reorder(sample_name, sample_number): sample_137<br />percent_aligned: 75.85013<br />as.factor(tank): 7","reorder(sample_name, sample_number): sample_138<br />percent_aligned: 77.93514<br />as.factor(tank): 7","reorder(sample_name, sample_number): sample_139<br />percent_aligned: 78.10972<br />as.factor(tank): 7","reorder(sample_name, sample_number): sample_140<br />percent_aligned: 75.89559<br />as.factor(tank): 7"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,190,103,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"7","legendgroup":"7","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0],"x":[75,76,77,78,79],"y":[74.739847670978776,71.385380896714764,58.313748092544493,69.230076755631458,79.324083185693055],"text":["reorder(sample_name, sample_number): sample_147<br />percent_aligned: 74.73985<br />as.factor(tank): 8","reorder(sample_name, sample_number): sample_148<br />percent_aligned: 71.38538<br />as.factor(tank): 8","reorder(sample_name, sample_number): sample_149<br />percent_aligned: 58.31375<br />as.factor(tank): 8","reorder(sample_name, sample_number): sample_150<br />percent_aligned: 69.23008<br />as.factor(tank): 8","reorder(sample_name, sample_number): sample_RESUB-156<br />percent_aligned: 79.32408<br />as.factor(tank): 8"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,193,154,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"8","legendgroup":"8","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0],"x":[40,41,42,43,44],"y":[75.229249487301502,72.720421898167999,74.399464671647593,73.692177330834483,75.472746916236488],"text":["reorder(sample_name, sample_number): sample_67<br />percent_aligned: 75.22925<br />as.factor(tank): 9","reorder(sample_name, sample_number): sample_68<br />percent_aligned: 72.72042<br />as.factor(tank): 9","reorder(sample_name, sample_number): sample_69<br />percent_aligned: 74.39946<br />as.factor(tank): 9","reorder(sample_name, sample_number): sample_70<br />percent_aligned: 73.69218<br />as.factor(tank): 9","reorder(sample_name, sample_number): sample_RESUB-76<br />percent_aligned: 75.47275<br />as.factor(tank): 9"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"9","legendgroup":"9","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000213,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0,0,0,0],"x":[32,33,34,35,36,37,38,39],"y":[78.373689421732564,75.29179932210657,39.498585386987685,65.78697556179219,74.564147342965597,68.48448055551313,79.012048738854773,73.069596027789444],"text":["reorder(sample_name, sample_number): sample_57<br />percent_aligned: 78.37369<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_57-G<br />percent_aligned: 75.29180<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_57-S<br />percent_aligned: 39.49859<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_58<br />percent_aligned: 65.78698<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_58-G<br />percent_aligned: 74.56415<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_58-S<br />percent_aligned: 68.48448<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_59<br />percent_aligned: 79.01205<br />as.factor(tank): 10","reorder(sample_name, sample_number): sample_60<br />percent_aligned: 73.06960<br />as.factor(tank): 10"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,184,231,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"10","legendgroup":"10","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0],"x":[28,29,30,31],"y":[76.052982059048475,75.920547971626178,72.836534408040066,74.506506231057386],"text":["reorder(sample_name, sample_number): sample_47<br />percent_aligned: 76.05298<br />as.factor(tank): 11","reorder(sample_name, sample_number): sample_48<br />percent_aligned: 75.92055<br />as.factor(tank): 11","reorder(sample_name, sample_number): sample_49<br />percent_aligned: 72.83653<br />as.factor(tank): 11","reorder(sample_name, sample_number): sample_50<br />percent_aligned: 74.50651<br />as.factor(tank): 11"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,169,255,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"11","legendgroup":"11","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0],"x":[23,24,25,26,27],"y":[78.137666740364068,72.289551451902852,75.65512701806442,70.292892735598585,63.887216917731735],"text":["reorder(sample_name, sample_number): sample_37<br />percent_aligned: 78.13767<br />as.factor(tank): 12","reorder(sample_name, sample_number): sample_38<br />percent_aligned: 72.28955<br />as.factor(tank): 12","reorder(sample_name, sample_number): sample_39<br />percent_aligned: 75.65513<br />as.factor(tank): 12","reorder(sample_name, sample_number): sample_40<br />percent_aligned: 70.29289<br />as.factor(tank): 12","reorder(sample_name, sample_number): sample_41<br />percent_aligned: 63.88722<br />as.factor(tank): 12"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(132,148,255,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"12","legendgroup":"12","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0],"x":[45,46,47,48],"y":[76.078769568721725,79.148045486106852,74.173493451808454,79.228327961793738],"text":["reorder(sample_name, sample_number): sample_78<br />percent_aligned: 76.07877<br />as.factor(tank): 13","reorder(sample_name, sample_number): sample_79<br />percent_aligned: 79.14805<br />as.factor(tank): 13","reorder(sample_name, sample_number): sample_80<br />percent_aligned: 74.17349<br />as.factor(tank): 13","reorder(sample_name, sample_number): sample_83<br />percent_aligned: 79.22833<br />as.factor(tank): 13"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(199,124,255,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"13","legendgroup":"13","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0],"x":[49,50,51,52],"y":[63.335167636042634,70.964840711590696,75.422840105451399,78.309635680699046],"text":["reorder(sample_name, sample_number): sample_88<br />percent_aligned: 63.33517<br />as.factor(tank): 14","reorder(sample_name, sample_number): sample_90<br />percent_aligned: 70.96484<br />as.factor(tank): 14","reorder(sample_name, sample_number): sample_91<br />percent_aligned: 75.42284<br />as.factor(tank): 14","reorder(sample_name, sample_number): sample_RESUB-94<br />percent_aligned: 78.30964<br />as.factor(tank): 14"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(237,104,237,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"14","legendgroup":"14","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0],"x":[53,54,55,56],"y":[78.01720630789174,69.647810861164388,71.94534440733608,79.27637489177765],"text":["reorder(sample_name, sample_number): sample_97<br />percent_aligned: 78.01721<br />as.factor(tank): 15","reorder(sample_name, sample_number): sample_98<br />percent_aligned: 69.64781<br />as.factor(tank): 15","reorder(sample_name, sample_number): sample_99<br />percent_aligned: 71.94534<br />as.factor(tank): 15","reorder(sample_name, sample_number): sample_100<br />percent_aligned: 79.27637<br />as.factor(tank): 15"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,97,204,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"15","legendgroup":"15","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568,0.90000000000000568],"base":[0,0,0,0,0],"x":[57,58,59,60,61],"y":[79.166188449495436,68.727363448190431,78.212098206040622,81.473563533017696,77.510652847664048],"text":["reorder(sample_name, sample_number): sample_107<br />percent_aligned: 79.16619<br />as.factor(tank): 16","reorder(sample_name, sample_number): sample_108<br />percent_aligned: 68.72736<br />as.factor(tank): 16","reorder(sample_name, sample_number): sample_109<br />percent_aligned: 78.21210<br />as.factor(tank): 16","reorder(sample_name, sample_number): sample_110<br />percent_aligned: 81.47356<br />as.factor(tank): 16","reorder(sample_name, sample_number): sample_RESUB-116<br />percent_aligned: 77.51065<br />as.factor(tank): 16"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,104,161,1)","line":{"width":1.8897637795275593,"color":"transparent"}},"name":"16","legendgroup":"16","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.228310502283104,"r":7.3059360730593621,"b":40.05685961550482,"l":37.260273972602747},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,79.599999999999994],"tickmode":"array","ticktext":["sample_1","sample_2","sample_3","sample_4","sample_5","sample_10","sample_11","sample_12","sample_13","sample_18","sample_19","sample_19-G","sample_19-S","sample_20","sample_20-G","sample_20-S","sample_21","sample_28","sample_29","sample_30","sample_31","sample_RESUB-36","sample_37","sample_38","sample_39","sample_40","sample_41","sample_47","sample_48","sample_49","sample_50","sample_57","sample_57-G","sample_57-S","sample_58","sample_58-G","sample_58-S","sample_59","sample_60","sample_67","sample_68","sample_69","sample_70","sample_RESUB-76","sample_78","sample_79","sample_80","sample_83","sample_88","sample_90","sample_91","sample_RESUB-94","sample_97","sample_98","sample_99","sample_100","sample_107","sample_108","sample_109","sample_110","sample_RESUB-116","sample_117","sample_118","sample_119","sample_120","sample_121","sample_127","sample_128","sample_129","sample_131","sample_137","sample_138","sample_139","sample_140","sample_147","sample_148","sample_149","sample_150","sample_RESUB-156"],"tickvals":[1,2,3,4,5,6,7,8,9,10,10.999999999999998,12,13,14,15,16,17,18,19,20,21,22,23,24.000000000000004,25,26,27,28,29,30,31.000000000000004,32,33,34,35,36,37,38,39,40,41.000000000000007,42,43,44,45,46.000000000000007,47,48,49,50,50.999999999999993,52,53,54,55.000000000000007,56,57,58,59,59.999999999999993,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79],"categoryorder":"array","categoryarray":["sample_1","sample_2","sample_3","sample_4","sample_5","sample_10","sample_11","sample_12","sample_13","sample_18","sample_19","sample_19-G","sample_19-S","sample_20","sample_20-G","sample_20-S","sample_21","sample_28","sample_29","sample_30","sample_31","sample_RESUB-36","sample_37","sample_38","sample_39","sample_40","sample_41","sample_47","sample_48","sample_49","sample_50","sample_57","sample_57-G","sample_57-S","sample_58","sample_58-G","sample_58-S","sample_59","sample_60","sample_67","sample_68","sample_69","sample_70","sample_RESUB-76","sample_78","sample_79","sample_80","sample_83","sample_88","sample_90","sample_91","sample_RESUB-94","sample_97","sample_98","sample_99","sample_100","sample_107","sample_108","sample_109","sample_110","sample_RESUB-116","sample_117","sample_118","sample_119","sample_120","sample_121","sample_127","sample_128","sample_129","sample_131","sample_137","sample_138","sample_139","sample_140","sample_147","sample_148","sample_149","sample_150","sample_RESUB-156"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":9.2984640929846396},"tickangle":-60,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"reorder(sample_name, sample_number)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-4.188944846742058,87.967841781583218],"tickmode":"array","ticktext":["0","20","40","60","80"],"tickvals":[0,20,40,60,80],"categoryorder":"array","categoryarray":["0","20","40","60","80"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"percent_aligned","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"title":{"text":"as.factor(tank)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"a1a87fffe88e":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"a1a87fffe88e","visdat":{"a1a87fffe88e":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

## 2.6 Trinity Matrix with Kallisto Output

``` bash
# Load bash variables into memory
source .bashvars

cd ${kallisto_output_dir}

${programs_array[trinity_abund_to_matrix]} \
--est_method 'kallisto' \
--gene_trans_map 'none' \
--out_prefix 'kallisto' \
--name_sample_by_basedir ${kallisto_output_dir}/kallisto_quant_*/abundance.tsv

ls -lh ${kallisto_output_dir}
```

    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_100/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_107/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_108/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_109/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_10/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_110/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_117/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_118/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_119/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_11/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_120/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_121/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_127/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_128/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_129/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_12/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_131/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_137/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_138/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_139/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_13/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_140/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_147/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_148/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_149/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_150/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_18/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19-G/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_19-S/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_1/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20-G/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_20-S/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_21/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_28/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_29/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_2/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_30/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_31/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_37/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_38/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_39/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_3/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_40/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_41/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_47/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_48/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_49/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_4/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_50/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57-G/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_57-S/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58-G/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_58-S/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_59/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_5/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_60/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_67/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_68/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_69/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_70/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_78/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_79/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_80/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_83/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_88/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_90/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_91/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_97/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_98/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_99/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-116/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-156/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-36/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-76/abundance.tsv
    -reading file: /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/06-cod-RNAseq-alignment/kallisto/kallisto_quant_RESUB-94/abundance.tsv


    * Outputting combined matrix.

    /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrixCMD: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2 
    sh: 1: R: not found
    Error, cmd: R --no-save --no-restore --no-site-file --no-init-file -q < kallisto.isoform.TPM.not_cross_norm.runTMM.R 1>&2  died with ret (32512)  at /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl line 105.
    Error, CMD: /home/shared/trinityrnaseq-v2.12.0/util/support_scripts/run_TMM_scale_matrix.pl --matrix kallisto.isoform.TPM.not_cross_norm > kallisto.isoform.TMM.EXPR.matrix died with ret 6400 at /home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl line 385.
    total 1.5G
    -rw-r--r-- 1 shedurkin labmembers 1.5G Mar 18 16:08 G_macrocephalus_kallisto_index.idx
    -rw-r--r-- 1 shedurkin labmembers  20M Apr 29 23:20 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Apr 29 23:20 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers  24M Apr 29 23:20 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Apr 29 23:20 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:46 kallisto_quant_1
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:11 kallisto_quant_10
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:03 kallisto_quant_100
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:04 kallisto_quant_100.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:05 kallisto_quant_107
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:05 kallisto_quant_107.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:07 kallisto_quant_108
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:08 kallisto_quant_108.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:09 kallisto_quant_109
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:10 kallisto_quant_109.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:12 kallisto_quant_10.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:22 kallisto_quant_11
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:13 kallisto_quant_110
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:14 kallisto_quant_110.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:15 kallisto_quant_117
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:16 kallisto_quant_117.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:18 kallisto_quant_118
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:18 kallisto_quant_118.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:20 kallisto_quant_119
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:20 kallisto_quant_119.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:22 kallisto_quant_11.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:33 kallisto_quant_12
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:24 kallisto_quant_120
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:24 kallisto_quant_120.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:26 kallisto_quant_121
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:27 kallisto_quant_121.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:28 kallisto_quant_127
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:29 kallisto_quant_127.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:30 kallisto_quant_128
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:30 kallisto_quant_128.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:31 kallisto_quant_129
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:31 kallisto_quant_129.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:33 kallisto_quant_12.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:39 kallisto_quant_13
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:34 kallisto_quant_131
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:34 kallisto_quant_131.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:35 kallisto_quant_137
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:35 kallisto_quant_137.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:36 kallisto_quant_138
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:36 kallisto_quant_138.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:37 kallisto_quant_139
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:37 kallisto_quant_139.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:39 kallisto_quant_13.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:40 kallisto_quant_140
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:40 kallisto_quant_140.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:41 kallisto_quant_147
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:41 kallisto_quant_147.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:43 kallisto_quant_148
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:44 kallisto_quant_148.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 15:51 kallisto_quant_149
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 15:52 kallisto_quant_149.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:52 kallisto_quant_150
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 22:53 kallisto_quant_150.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:37 kallisto_quant_18
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:37 kallisto_quant_18.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:41 kallisto_quant_19
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:39 kallisto_quant_19-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:40 kallisto_quant_19-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:41 kallisto_quant_19.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:43 kallisto_quant_19-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:44 kallisto_quant_19-S.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:46 kallisto_quant_1.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:59 kallisto_quant_2
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:50 kallisto_quant_20
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:48 kallisto_quant_20-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:49 kallisto_quant_20-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:50 kallisto_quant_20.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:52 kallisto_quant_20-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:52 kallisto_quant_20-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:54 kallisto_quant_21
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:54 kallisto_quant_21.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:56 kallisto_quant_28
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:56 kallisto_quant_28.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 17:57 kallisto_quant_29
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 17:58 kallisto_quant_29.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:00 kallisto_quant_2.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:53 kallisto_quant_3
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 18:01 kallisto_quant_30
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:02 kallisto_quant_30.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 18:03 kallisto_quant_31
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 18:04 kallisto_quant_31.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:47 kallisto_quant_37
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:48 kallisto_quant_37.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:49 kallisto_quant_38
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:49 kallisto_quant_38.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:51 kallisto_quant_39
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:51 kallisto_quant_39.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:53 kallisto_quant_3.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:04 kallisto_quant_4
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:54 kallisto_quant_40
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:55 kallisto_quant_40.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:56 kallisto_quant_41
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:56 kallisto_quant_41.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 20:58 kallisto_quant_47
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 20:58 kallisto_quant_47.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:00 kallisto_quant_48
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:00 kallisto_quant_48.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:02 kallisto_quant_49
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:02 kallisto_quant_49.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:04 kallisto_quant_4.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:26 kallisto_quant_5
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:06 kallisto_quant_50
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:06 kallisto_quant_50.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:14 kallisto_quant_57
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:10 kallisto_quant_57-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:11 kallisto_quant_57-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:14 kallisto_quant_57.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:15 kallisto_quant_57-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:15 kallisto_quant_57-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:20 kallisto_quant_58
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:18 kallisto_quant_58-G
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:18 kallisto_quant_58-G.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:20 kallisto_quant_58.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:22 kallisto_quant_58-S
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:22 kallisto_quant_58-S.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:24 kallisto_quant_59
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:24 kallisto_quant_59.log
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:26 kallisto_quant_5.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:28 kallisto_quant_60
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:28 kallisto_quant_60.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:30 kallisto_quant_67
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:30 kallisto_quant_67.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:32 kallisto_quant_68
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:32 kallisto_quant_68.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:34 kallisto_quant_69
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:34 kallisto_quant_69.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:35 kallisto_quant_70
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:36 kallisto_quant_70.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:37 kallisto_quant_78
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:37 kallisto_quant_78.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:39 kallisto_quant_79
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:39 kallisto_quant_79.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:40 kallisto_quant_80
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:41 kallisto_quant_80.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:42 kallisto_quant_83
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:43 kallisto_quant_83.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:44 kallisto_quant_88
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:44 kallisto_quant_88.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:46 kallisto_quant_90
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:46 kallisto_quant_90.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:48 kallisto_quant_91
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:48 kallisto_quant_91.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:50 kallisto_quant_97
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:50 kallisto_quant_97.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:52 kallisto_quant_98
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:52 kallisto_quant_98.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:53 kallisto_quant_99
    -rw-r--r-- 1 shedurkin labmembers 5.2K Apr 29 21:54 kallisto_quant_99.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:55 kallisto_quant_RESUB-116
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:55 kallisto_quant_RESUB-116.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:57 kallisto_quant_RESUB-156
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:57 kallisto_quant_RESUB-156.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 21:58 kallisto_quant_RESUB-36
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 21:59 kallisto_quant_RESUB-36.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:54 kallisto_quant_RESUB-76
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 22:54 kallisto_quant_RESUB-76.log
    drwxr-xr-x 2 shedurkin labmembers 4.0K Apr 29 22:13 kallisto_quant_RESUB-94
    -rw-r--r-- 1 shedurkin labmembers 5.3K Apr 29 22:14 kallisto_quant_RESUB-94.log
