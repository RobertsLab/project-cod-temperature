05-cod-RNAseq-trimming
================
Kathleen Durkin
2024-03-04

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-download-raw-rnaseq-reads"
  id="toc-2-download-raw-rnaseq-reads">2 Download raw RNAseq reads</a>
  - <a href="#21-verify-raw-read-checksums"
    id="toc-21-verify-raw-read-checksums">2.1 Verify raw read checksums</a>
- <a href="#3-fastqcmultiqc-on-raw-reads"
  id="toc-3-fastqcmultiqc-on-raw-reads">3 FastQC/MultiQC on raw reads</a>
- <a href="#4-trimming-with-flexbar" id="toc-4-trimming-with-flexbar">4
  Trimming with flexbar</a>
- <a href="#5-fastqcmultiqc-on-trimmed-reads"
  id="toc-5-fastqcmultiqc-on-trimmed-reads">5 FastQC/MultiQC on trimmed
  reads</a>
- <a href="#6-summary" id="toc-6-summary">6 Summary</a>

Code for trimming and QCing RNAseq data, to be used on [Pacific cod
RNAseq
data](https://shedurkin.github.io/Roberts-LabNotebook/posts/projects/pacific_cod/2023_12_13_pacific_cod.html),
raw reads found
[here](https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/)

FastQC/MultiQC assessment of raw and
[flexbar](https://github.com/seqan/flexbar)-trimmed sequences of

Inputs:

- RNAseq gzipped FastQs (e.g. `*.fastq.gz`)

Outputs:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  HTML reports for raw and trimmed reads.

- [`MultiQC`](https://multiqc.info/) HTML summaries of
  [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  for raw and trimmed reads.

- Trimmed reads: `*flexbar_trim.fastq.gz`

(info about library prep/sequencing)

------------------------------------------------------------------------

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export cod_dir=/home/shared/8TB_HDD_02/shedurkin/project-cod-temperature'
echo 'export output_dir_top=${cod_dir}/output/05-cod-RNAseq-trimming'
echo 'export raw_fastqc_dir=${output_dir_top}/raw-fastqc'
echo 'export raw_reads_dir=${cod_dir}/data/05-cod-RNAseq-trimming/raw-reads'
echo 'export raw_reads_url="https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/"'
echo 'export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc'
echo 'export trimmed_reads_dir=${output_dir_top}/trimmed-reads'
echo ""

echo "# Paths to programs"
echo 'export fastqc=/home/shared/FastQC-0.12.1/fastqc'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo 'export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar'
echo ""

echo "# Set FastQ filename patterns"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1_*.fastq.gz'"
echo "export R2_fastq_pattern='*_R2_*.fastq.gz'"
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""

echo "# Input/output files"
echo 'export fastq_checksums=input_fastq_checksums.md5'
echo 'export trimmed_checksums=trimmed_fastq_checksums.md5'
echo 'export NEB_adapters_fasta=NEB-adapters.fasta'
echo ""

echo "## Illumina adapters"
echo 'export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"'
echo 'export second_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"'
echo ""

echo "## Inititalize arrays"
echo 'export fastq_array_R1=()'
echo 'export fastq_array_R2=()'
echo 'export raw_fastqs_array=()'
echo 'export R1_names_array=()'
echo 'export R2_names_array=()'
echo 'export trimmed_fastqs_array=()'
echo ""

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[fastqc]="${fastqc}" \'
echo '[multiqc]="${multiqc}" \'
echo '[flexbar]="${flexbar}"'
echo ")"
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export cod_dir=/home/shared/8TB_HDD_02/shedurkin/project-cod-temperature
    export output_dir_top=${cod_dir}/output/05-cod-RNAseq-trimming
    export raw_fastqc_dir=${output_dir_top}/raw-fastqc
    export raw_reads_dir=${cod_dir}/data/05-cod-RNAseq-trimming/raw-reads
    export raw_reads_url="https://owl.fish.washington.edu/nightingales/G_macrocephalus/30-943133806/"
    export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc
    export trimmed_reads_dir=${output_dir_top}/trimmed-reads

    # Paths to programs
    export fastqc=/home/shared/FastQC-0.12.1/fastqc
    export multiqc=/home/sam/programs/mambaforge/bin/multiqc
    export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar

    # Set FastQ filename patterns
    export fastq_pattern='*.fastq.gz'
    export R1_fastq_pattern='*_R1_*.fastq.gz'
    export R2_fastq_pattern='*_R2_*.fastq.gz'

    # Set number of CPUs to use
    export threads=20

    # Input/output files
    export fastq_checksums=input_fastq_checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5
    export NEB_adapters_fasta=NEB-adapters.fasta

    ## Illumina adapters
    export first_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    export second_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

    ## Inititalize arrays
    export fastq_array_R1=()
    export fastq_array_R2=()
    export raw_fastqs_array=()
    export R1_names_array=()
    export R2_names_array=()
    export trimmed_fastqs_array=()

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [fastqc]="${fastqc}" \
    [multiqc]="${multiqc}" \
    [flexbar]="${flexbar}"
    )

# 2 Download raw RNAseq reads

Reads are downloaded from:

The `--cut-dirs 3` command cuts the preceding directory structure
(i.e. `nightingales/P_meandrina/30-852430235/`) so that we just end up
with the reads.

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--accept ${fastq_pattern} ${raw_reads_url}

ls -lh "${raw_reads_dir}"
```

## 2.1 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

cd "${raw_reads_dir}"

md5sum checksums.md5 --check
```

# 3 FastQC/MultiQC on raw reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############


# Create array of raw FastQs
raw_fastqs_array=(${raw_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
raw_fastqc_list=$(echo "${raw_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
### NOTE: Do NOT quote raw_fastqc_list
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${raw_fastqc_dir} \
--quiet \
${raw_fastqc_list}

echo "FastQC on raw reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${raw_fastqc_dir} -o ${raw_fastqc_dir}

echo ""
echo "MultiQC on raw FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${raw_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""
```

``` bash
# View directory contents
ls -lh ${raw_fastqc_dir}
```

    total 72M
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:06 01_temp-size-analysis_files
    -rw-r--r-- 1 shedurkin labmembers 1.4M May 23 15:06 01_temp-size-analysis.html
    -rw-r--r-- 1 shedurkin labmembers  17K May 23 15:06 01_temp-size-analysis.md
    -rw-r--r-- 1 shedurkin labmembers 361K Oct 25  2023 01_temp-size-analysis.pdf
    -rw-r--r-- 1 shedurkin labmembers 9.7K May 23 15:06 01_temp-size-analysis.Rmd
    drwxr-xr-x 4 shedurkin labmembers 4.0K May  8 11:05 03.2-genome-annotation_cache
    drwxr-xr-x 4 shedurkin labmembers 4.0K May  8 11:07 03.2-genome-annotation_files
    -rw-r--r-- 1 shedurkin labmembers 1.1M May 22 14:45 03.2-genome-annotation.html
    -rw-r--r-- 1 shedurkin labmembers 5.0K May 22 14:45 03.2-genome-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 4.2K May 22 14:45 03.2-genome-annotation.Rmd
    drwxr-xr-x 4 shedurkin labmembers 4.0K Apr 17 13:00 03-transcriptome-annotation_cache
    drwxr-xr-x 4 shedurkin labmembers 4.0K Apr 17 13:00 03-transcriptome-annotation_files
    -rw-r--r-- 1 shedurkin labmembers  50M Apr 24 12:59 03-transcriptome-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 9.6K Apr 24 12:54 03-transcriptome-annotation.Rmd
    -rw-r--r-- 1 shedurkin labmembers 612K Nov  2  2023 04-RNASeq-sample-size.html
    -rw-r--r-- 1 shedurkin labmembers 3.0K Nov  2  2023 04-RNASeq-sample-size.md
    -rw-r--r-- 1 shedurkin labmembers 2.5K Nov  2  2023 04-RNASeq-sample-size.Rmd
    -rw-r--r-- 1 shedurkin labmembers 1.1M Mar 22 12:18 05-cod-RNAseq-trimming.html
    -rw-r--r-- 1 shedurkin labmembers  15K Mar 22 12:15 05-cod-RNAseq-trimming.md
    -rw-r--r-- 1 shedurkin labmembers  11K Jun 10 12:20 05-cod-RNAseq-trimming.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  3 14:59 06.2-cod-RNAseq-alignment-genome_files
    -rw-r--r-- 1 shedurkin labmembers 1.8M May  3 14:59 06.2-cod-RNAseq-alignment-genome.html
    -rw-r--r-- 1 shedurkin labmembers  36K May  3 14:59 06.2-cod-RNAseq-alignment-genome.md
    -rw-r--r-- 1 shedurkin labmembers  16K May 21 14:01 06.2-cod-RNAseq-alignment-genome.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  3 15:02 06-cod-RNAseq-alignment_files
    -rw-r--r-- 1 shedurkin labmembers 1.6M May  3 15:02 06-cod-RNAseq-alignment.html
    -rw-r--r-- 1 shedurkin labmembers  52K May  3 15:01 06-cod-RNAseq-alignment.md
    -rw-r--r-- 1 shedurkin labmembers  11K May  3 15:01 06-cod-RNAseq-alignment.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  8 13:42 07.2.1-cod-RNAseq-DESeq2-genome-exon_files
    -rw-r--r-- 1 shedurkin labmembers  71K May  8 13:48 07.2.1-cod-RNAseq-DESeq2-genome-exon.md
    -rw-r--r-- 1 shedurkin labmembers  33K May  8 13:43 07.2.1-cod-RNAseq-DESeq2-genome-exon.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  8 13:28 07.2.2-cod-RNAseq-DESeq2-genome-gene_files
    -rw-r--r-- 1 shedurkin labmembers 4.6M May  8 13:28 07.2.2-cod-RNAseq-DESeq2-genome-gene.html
    -rw-r--r-- 1 shedurkin labmembers  71K May  8 13:23 07.2.2-cod-RNAseq-DESeq2-genome-gene.md
    -rw-r--r-- 1 shedurkin labmembers  33K May  8 13:23 07.2.2-cod-RNAseq-DESeq2-genome-gene.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 21 15:24 07-cod-RNAseq-DESeq2_files
    -rw-r--r-- 1 shedurkin labmembers 4.7M May 21 15:24 07-cod-RNAseq-DESeq2.html
    -rw-r--r-- 1 shedurkin labmembers  71K May 21 15:16 07-cod-RNAseq-DESeq2.md
    -rw-r--r-- 1 shedurkin labmembers  35K Jun 10 12:16 07-cod-RNAseq-DESeq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon_files
    -rw-r--r-- 1 shedurkin labmembers 1.4M May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.html
    -rw-r--r-- 1 shedurkin labmembers 159K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.md
    -rw-r--r-- 1 shedurkin labmembers  17K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 13 13:29 08-cod-RNAseq-GO-annotation_files
    -rw-r--r-- 1 shedurkin labmembers 1.3M May 13 13:29 08-cod-RNAseq-GO-annotation.html
    -rw-r--r-- 1 shedurkin labmembers 129K May 13 13:28 08-cod-RNAseq-GO-annotation.md
    -rw-r--r-- 1 shedurkin labmembers  18K May 13 13:29 08-cod-RNAseq-GO-annotation.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 22 11:36 09-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers 8.4K May 22 11:36 09-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 1.8K May 22 11:36 09-hisat-deseq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10.1-hisat-deseq2_cache
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10.1-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers  19K May 23 15:08 10.1-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 7.6K May 23 15:08 10.1-hisat-deseq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10-hisat-deseq2_cache
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers  18K May 23 15:08 10-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 7.8K May 23 15:08 10-hisat-deseq2.Rmd
    -rw-r--r-- 1 shedurkin labmembers 112K May 23 15:08 11-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 4.4K May 23 15:08 11-annotation.Rmd
    -rw-r--r-- 1 shedurkin labmembers 2.1K May 23 15:08 12-DEG-enrich.Rmd
    -rw-r--r-- 1 shedurkin labmembers  844 Jun  6 14:42 13.0.0-RNAseq-edgeR.Rmd
    -rw-r--r-- 1 shedurkin labmembers 1.6K May 22 11:36 README.md
    drwxr-xr-x 3 shedurkin labmembers 4.0K Oct 25  2023 rsconnect
    -rw-r--r-- 1 shedurkin labmembers 1.7M Oct 25  2023 Temp-Size-Graph-Box-Line.html
    -rw-r--r-- 1 shedurkin labmembers 3.7K Oct 25  2023 Temp-Size-Graph-Box-Line.log
    -rw-r--r-- 1 shedurkin labmembers 2.7K Oct 25  2023 Temp-Size-Graph-Box-Line.Rmd
    -rw-r--r-- 1 shedurkin labmembers  14K Oct 25  2023 Temp-Size-Graph-Box-Line.tex

Samples 10, 13, 19, 37, 41, 48, 98, 129, 149, 57-S, 58-S had low
quantity/quality of RNA after RNA extraction (see [Azenta RNA extraction
QC
report](https://github.com/RobertsLab/project-cod-temperature/blob/main/data/Sample.QC.report.of_30-943133806_240118025106.pdf))

**Notes:** - Not enough starting material: 19, 41, 48, 98, 129, 149, -
RIN \< 6: 10, 13, 37, 57-S, 58-S - DV200 \< 70: 19, 41, 48, 98, 129,
149, 57-S

RNA Integrity Number (RIN) is a measure of RNA quality ranging from 1 -
10, with 1 representing highly degraded RNA and 10 representing intact
RNA. Most researchers aim for RIN values of \> 8, but Azenta considers
anything \> 6 acceptable. It is also generally assumed that RIN is
representative of mRNA quality. It’s worth noting there is some debate
surrounding the accuracy of RIN as a measure of RNA integrity – see
[ScienceDirect’s summary of
RIN](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/rna-integrity-number)

DV200 is another measure of RNA quality that represents the percentage
of RNA fragments of \>200 nucleotides in length. A recent study
suggested that DV200 may be a more meaningful assessment of mRNA quality
than RIN ([Matsubara et al 2020](https://doi.org/10.1155/2020/9349132))
for NGS, but both are still commonly used.

# 4 Trimming with [flexbar](https://github.com/seqan/flexbar)

``` bash
# Load bash variables into memory
source .bashvars

# Change to directory with raw reads
cd "${raw_reads_dir}"

# Create arrays of FastQ R1 files and sample names
# Do NOT quote R1_fastq_pattern variable
for fastq in ${R1_fastq_pattern}
do
  fastq_array_R1+=("${fastq}")

  # Use parameter substitution to remove all text up to and including last "." from
  # right side of string.
  R1_names_array+=("${fastq%%.*}")
done

# Create array of FastQ R2 files
# Do NOT quote R2_fastq_pattern variable
for fastq in ${R2_fastq_pattern}
do
  fastq_array_R2+=("${fastq}")

  # Use parameter substitution to remove all text up to and including last "." from
  # right side of string.
  R2_names_array+=("${fastq%%.*}")
done


############ RUN FLEXBAR ############
# Uses parameter substitution (e.g. ${R1_sample_name%%_*})to rm the _R[12]
# Uses built-in adapter presets (https://github.com/seqan/flexbar/wiki/Manual#adapter-removal-presets)
# --adapter-preset TruSeq: several adapter presets for Illumina libraries are included in Flexbar
# --adapter-pair-overlap ON: Recommended by NEB sRNA kit
# --qtrim-threshold 25: Minimum quality
# --qtrim-format i1.8: Sets sequencer as Illumina
# --target: Sets file naming patterns
# --zip-output GZ: Sets type of compression. GZ = gzip
#

# Run flexbar on files
echo "Beginning flexbar trimming."
echo ""

time \
for index in "${!fastq_array_R1[@]}"
do
  R1_sample_name="${R1_names_array[index]}"
  R2_sample_name="${R2_names_array[index]}"

  # Begin flexbar trimming
  ${programs_array[flexbar]} \
  --reads ${fastq_array_R1[index]} \
  --reads2 ${fastq_array_R2[index]}  \
  --adapter-preset "TruSeq" \
  --adapter-pair-overlap ON \
  --qtrim-format i1.8 \
  --qtrim-threshold 25 \
  --threads ${threads} \
  --target "${trimmed_reads_dir}/${R1_sample_name%%_*}.flexbar_trim.R" \
  --zip-output GZ
        
    # Move to trimmed directory
    # This is done so checksums file doesn't include excess path
    cd ${trimmed_reads_dir}

    # Generate md5 checksums for newly trimmed files
    {
      md5sum "${R1_sample_name%%_*}.flexbar_trim.R_1.fastq.gz"
      md5sum "${R2_sample_name%%_*}.flexbar_trim.R_2.fastq.gz"
    } >> "${trimmed_checksums}"
    
    # Change back to to raw reads directory
    cd "${raw_reads_dir}"

done

echo ""
echo "flexbar trimming complete."
echo ""

echo "Trimmed FastQs MD5 checksums:"
echo ""

cat "${trimmed_reads_dir}/${trimmed_checksums}"

############ END FLEXBAR ############
```

# 5 FastQC/MultiQC on trimmed reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############

### NOTE: Do NOT quote raw_fastqc_list
# Create array of trimmed FastQs
trimmed_fastqs_array=(${trimmed_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
trimmed_fastqc_list=$(echo "${trimmed_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${trimmed_fastqc_dir} \
--quiet \
${trimmed_fastqc_list}

echo "FastQC on trimmed reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${trimmed_fastqc_dir} -o ${trimmed_fastqc_dir}

echo ""
echo "MultiQC on trimmed FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${trimmed_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""
```

``` bash
# View directory contents
ls -lh ${trimmed_fastqc_dir}
```

    total 72M
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:06 01_temp-size-analysis_files
    -rw-r--r-- 1 shedurkin labmembers 1.4M May 23 15:06 01_temp-size-analysis.html
    -rw-r--r-- 1 shedurkin labmembers  17K May 23 15:06 01_temp-size-analysis.md
    -rw-r--r-- 1 shedurkin labmembers 361K Oct 25  2023 01_temp-size-analysis.pdf
    -rw-r--r-- 1 shedurkin labmembers 9.7K May 23 15:06 01_temp-size-analysis.Rmd
    drwxr-xr-x 4 shedurkin labmembers 4.0K May  8 11:05 03.2-genome-annotation_cache
    drwxr-xr-x 4 shedurkin labmembers 4.0K May  8 11:07 03.2-genome-annotation_files
    -rw-r--r-- 1 shedurkin labmembers 1.1M May 22 14:45 03.2-genome-annotation.html
    -rw-r--r-- 1 shedurkin labmembers 5.0K May 22 14:45 03.2-genome-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 4.2K May 22 14:45 03.2-genome-annotation.Rmd
    drwxr-xr-x 4 shedurkin labmembers 4.0K Apr 17 13:00 03-transcriptome-annotation_cache
    drwxr-xr-x 4 shedurkin labmembers 4.0K Apr 17 13:00 03-transcriptome-annotation_files
    -rw-r--r-- 1 shedurkin labmembers  50M Apr 24 12:59 03-transcriptome-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 9.6K Apr 24 12:54 03-transcriptome-annotation.Rmd
    -rw-r--r-- 1 shedurkin labmembers 612K Nov  2  2023 04-RNASeq-sample-size.html
    -rw-r--r-- 1 shedurkin labmembers 3.0K Nov  2  2023 04-RNASeq-sample-size.md
    -rw-r--r-- 1 shedurkin labmembers 2.5K Nov  2  2023 04-RNASeq-sample-size.Rmd
    -rw-r--r-- 1 shedurkin labmembers 1.1M Mar 22 12:18 05-cod-RNAseq-trimming.html
    -rw-r--r-- 1 shedurkin labmembers  15K Mar 22 12:15 05-cod-RNAseq-trimming.md
    -rw-r--r-- 1 shedurkin labmembers  11K Jun 10 12:20 05-cod-RNAseq-trimming.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  3 14:59 06.2-cod-RNAseq-alignment-genome_files
    -rw-r--r-- 1 shedurkin labmembers 1.8M May  3 14:59 06.2-cod-RNAseq-alignment-genome.html
    -rw-r--r-- 1 shedurkin labmembers  36K May  3 14:59 06.2-cod-RNAseq-alignment-genome.md
    -rw-r--r-- 1 shedurkin labmembers  16K May 21 14:01 06.2-cod-RNAseq-alignment-genome.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  3 15:02 06-cod-RNAseq-alignment_files
    -rw-r--r-- 1 shedurkin labmembers 1.6M May  3 15:02 06-cod-RNAseq-alignment.html
    -rw-r--r-- 1 shedurkin labmembers  52K May  3 15:01 06-cod-RNAseq-alignment.md
    -rw-r--r-- 1 shedurkin labmembers  11K May  3 15:01 06-cod-RNAseq-alignment.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  8 13:42 07.2.1-cod-RNAseq-DESeq2-genome-exon_files
    -rw-r--r-- 1 shedurkin labmembers  71K May  8 13:48 07.2.1-cod-RNAseq-DESeq2-genome-exon.md
    -rw-r--r-- 1 shedurkin labmembers  33K May  8 13:43 07.2.1-cod-RNAseq-DESeq2-genome-exon.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May  8 13:28 07.2.2-cod-RNAseq-DESeq2-genome-gene_files
    -rw-r--r-- 1 shedurkin labmembers 4.6M May  8 13:28 07.2.2-cod-RNAseq-DESeq2-genome-gene.html
    -rw-r--r-- 1 shedurkin labmembers  71K May  8 13:23 07.2.2-cod-RNAseq-DESeq2-genome-gene.md
    -rw-r--r-- 1 shedurkin labmembers  33K May  8 13:23 07.2.2-cod-RNAseq-DESeq2-genome-gene.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 21 15:24 07-cod-RNAseq-DESeq2_files
    -rw-r--r-- 1 shedurkin labmembers 4.7M May 21 15:24 07-cod-RNAseq-DESeq2.html
    -rw-r--r-- 1 shedurkin labmembers  71K May 21 15:16 07-cod-RNAseq-DESeq2.md
    -rw-r--r-- 1 shedurkin labmembers  35K Jun 10 12:16 07-cod-RNAseq-DESeq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon_files
    -rw-r--r-- 1 shedurkin labmembers 1.4M May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.html
    -rw-r--r-- 1 shedurkin labmembers 159K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.md
    -rw-r--r-- 1 shedurkin labmembers  17K May 22 14:44 08.2.1-cod-RNAseq-GO-annotation-genome-exon.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 13 13:29 08-cod-RNAseq-GO-annotation_files
    -rw-r--r-- 1 shedurkin labmembers 1.3M May 13 13:29 08-cod-RNAseq-GO-annotation.html
    -rw-r--r-- 1 shedurkin labmembers 129K May 13 13:28 08-cod-RNAseq-GO-annotation.md
    -rw-r--r-- 1 shedurkin labmembers  18K May 13 13:29 08-cod-RNAseq-GO-annotation.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 22 11:36 09-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers 8.4K May 22 11:36 09-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 1.8K May 22 11:36 09-hisat-deseq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10.1-hisat-deseq2_cache
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10.1-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers  19K May 23 15:08 10.1-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 7.6K May 23 15:08 10.1-hisat-deseq2.Rmd
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10-hisat-deseq2_cache
    drwxr-xr-x 3 shedurkin labmembers 4.0K May 23 15:08 10-hisat-deseq2_files
    -rw-r--r-- 1 shedurkin labmembers  18K May 23 15:08 10-hisat-deseq2.md
    -rw-r--r-- 1 shedurkin labmembers 7.8K May 23 15:08 10-hisat-deseq2.Rmd
    -rw-r--r-- 1 shedurkin labmembers 112K May 23 15:08 11-annotation.md
    -rw-r--r-- 1 shedurkin labmembers 4.4K May 23 15:08 11-annotation.Rmd
    -rw-r--r-- 1 shedurkin labmembers 2.1K May 23 15:08 12-DEG-enrich.Rmd
    -rw-r--r-- 1 shedurkin labmembers  844 Jun  6 14:42 13.0.0-RNAseq-edgeR.Rmd
    -rw-r--r-- 1 shedurkin labmembers 1.6K May 22 11:36 README.md
    drwxr-xr-x 3 shedurkin labmembers 4.0K Oct 25  2023 rsconnect
    -rw-r--r-- 1 shedurkin labmembers 1.7M Oct 25  2023 Temp-Size-Graph-Box-Line.html
    -rw-r--r-- 1 shedurkin labmembers 3.7K Oct 25  2023 Temp-Size-Graph-Box-Line.log
    -rw-r--r-- 1 shedurkin labmembers 2.7K Oct 25  2023 Temp-Size-Graph-Box-Line.Rmd
    -rw-r--r-- 1 shedurkin labmembers  14K Oct 25  2023 Temp-Size-Graph-Box-Line.tex

# 6 Summary
