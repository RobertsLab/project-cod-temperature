06-cod-RNAseq-alignment-template
================
Kathleen Durkin
2024-01-12

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

# Create a Bash variables file

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
echo 'export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl'
echo ""


echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""


echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[kallisto]="${kallisto}" \'
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
    export trinity_abund_to_matrix=/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl

    # Set number of CPUs to use
    export threads=20

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [kallisto]="${kallisto}" \
    [trinity_abund_to_matrix]="${trinity_abund_to_matrix}" \
    )

# Align to reference transcriptome (Kallisto pseudoalignment)

## Retrieving the reference transcriptome

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

    total 565M
    drwxr-xr-x 3 shedurkin labmembers 4.0K Mar  4 11:05 05-cod-RNAseq-trimming
    -rw-r--r-- 1 shedurkin labmembers  13K Dec 27 15:45 Cod_RNAseq_NGS_Template_File.xlsx
    -rw-r--r-- 1 shedurkin labmembers  38M Oct 25 13:07 Gadus_macrocephalus.coding.gene.V1.cds
    -rw-r--r-- 1 shedurkin labmembers 169M Oct 16 19:28 GCF_031168955.1_ASM3116895v1_rna.fna
    -rw-r--r-- 1 shedurkin labmembers  47K Oct 25 12:54 Pcod Temp Growth experiment 2022-23 DATA.xlsx
    -rw-r--r-- 1 shedurkin labmembers 231K Mar  4 17:41 Sample.QC.report.of_30-943133806_240118025106.pdf
    -rw-r--r-- 1 shedurkin labmembers  12K Mar  4 17:41 Sample.QC.report.of_30-943133806_240118025106.xlsx
    -rw-r--r-- 1 shedurkin labmembers  12K Oct 25 12:54 temp-experiment.csv
    -rw-r--r-- 1 shedurkin labmembers 271M Oct 25 13:13 uniprot_sprot_r2023_04.fasta
    -rw-r--r-- 1 shedurkin labmembers  88M Oct 25 13:13 uniprot_sprot_r2023_04.fasta.gz

## Verify transcriptome FastA MD5 checksum

``` bash
# Load bash variables into memory
source .bashvars

cd "${transcriptome_fasta_dir}"

# Checksums file contains other files, so this just looks for the sRNAseq files.
md5sum --check <<< "${transcriptome_checksum}  ${transcriptome_fasta_name}.fna"
```

    GCF_031168955.1_ASM3116895v1_rna.fna: OK

## Building Index

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
    -rw-r--r-- 1 shedurkin labmembers  20M Mar 19 14:33 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Mar 19 14:33 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers  24M Mar 19 14:33 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Mar 19 14:33 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:00 kallisto_quant_1
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:32 kallisto_quant_10
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:28 kallisto_quant_100
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:29 kallisto_quant_107
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:30 kallisto_quant_108
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:31 kallisto_quant_109
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:37 kallisto_quant_11
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:33 kallisto_quant_110
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:34 kallisto_quant_117
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:35 kallisto_quant_118
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:37 kallisto_quant_119
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:44 kallisto_quant_12
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:39 kallisto_quant_120
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:39 kallisto_quant_121
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:41 kallisto_quant_127
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:41 kallisto_quant_128
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:42 kallisto_quant_129
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:47 kallisto_quant_13
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:43 kallisto_quant_131
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:44 kallisto_quant_137
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:45 kallisto_quant_138
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:45 kallisto_quant_139
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:47 kallisto_quant_140
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:48 kallisto_quant_147
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:49 kallisto_quant_148
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:56 kallisto_quant_149
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:52 kallisto_quant_150
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:54 kallisto_quant_18
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:58 kallisto_quant_19
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:57 kallisto_quant_19-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:59 kallisto_quant_19-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:09 kallisto_quant_2
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:02 kallisto_quant_20
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:02 kallisto_quant_20-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:04 kallisto_quant_20-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:05 kallisto_quant_21
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:07 kallisto_quant_28
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:07 kallisto_quant_29
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:15 kallisto_quant_3
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:09 kallisto_quant_30
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:11 kallisto_quant_31
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:11 kallisto_quant_37
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:13 kallisto_quant_38
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:13 kallisto_quant_39
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:21 kallisto_quant_4
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:15 kallisto_quant_40
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:17 kallisto_quant_41
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:17 kallisto_quant_47
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:19 kallisto_quant_48
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:19 kallisto_quant_49
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:32 kallisto_quant_5
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:22 kallisto_quant_50
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:24 kallisto_quant_57
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:24 kallisto_quant_57-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:25 kallisto_quant_57-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:27 kallisto_quant_58
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:27 kallisto_quant_58-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:30 kallisto_quant_58-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:30 kallisto_quant_59
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:32 kallisto_quant_60
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:34 kallisto_quant_67
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:35 kallisto_quant_68
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:36 kallisto_quant_69
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:36 kallisto_quant_70
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:38 kallisto_quant_78
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:38 kallisto_quant_79
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:40 kallisto_quant_80
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:40 kallisto_quant_83
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:42 kallisto_quant_88
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:42 kallisto_quant_90
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:44 kallisto_quant_91
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:44 kallisto_quant_97
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:46 kallisto_quant_98
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:46 kallisto_quant_99
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:48 kallisto_quant_RESUB-116
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:48 kallisto_quant_RESUB-156
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:49 kallisto_quant_RESUB-36
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:50 kallisto_quant_RESUB-76
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:51 kallisto_quant_RESUB-94

## Sample Quantification

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
        ${raw_reads_dir}/${R1_fastq} ${raw_reads_dir}/${R2_fastq}
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
have (#rawreads/2)+1 – one output folder for each pair of reads plus the
single kallisto index file

``` bash
# Load bash variables into memory
source .bashvars

# Count number of raw read files
cd ${raw_reads_dir}
echo "Number of raw reads:"
ls -1 | wc -l

# Count number of kallisto output 
cd ${kallisto_output_dir}
echo "Number of output files/folders"
ls -1 | wc -l
```

    Number of raw reads:
    158
    Number of output files/folders
    84

## Trinity Matrix with Kallisto Output

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
    -rw-r--r-- 1 shedurkin labmembers  20M Mar 19 14:43 kallisto.isoform.counts.matrix
    -rw-r--r-- 1 shedurkin labmembers    0 Mar 19 14:43 kallisto.isoform.TMM.EXPR.matrix
    -rw-r--r-- 1 shedurkin labmembers  24M Mar 19 14:43 kallisto.isoform.TPM.not_cross_norm
    -rw-r--r-- 1 shedurkin labmembers  532 Mar 19 14:43 kallisto.isoform.TPM.not_cross_norm.runTMM.R
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:00 kallisto_quant_1
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:32 kallisto_quant_10
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:28 kallisto_quant_100
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:29 kallisto_quant_107
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:30 kallisto_quant_108
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:31 kallisto_quant_109
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:37 kallisto_quant_11
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:33 kallisto_quant_110
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:34 kallisto_quant_117
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:35 kallisto_quant_118
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:37 kallisto_quant_119
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:44 kallisto_quant_12
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:39 kallisto_quant_120
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:39 kallisto_quant_121
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:41 kallisto_quant_127
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:41 kallisto_quant_128
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:42 kallisto_quant_129
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:47 kallisto_quant_13
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:43 kallisto_quant_131
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:44 kallisto_quant_137
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:45 kallisto_quant_138
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:45 kallisto_quant_139
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:47 kallisto_quant_140
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:48 kallisto_quant_147
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:49 kallisto_quant_148
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:56 kallisto_quant_149
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:52 kallisto_quant_150
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:54 kallisto_quant_18
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:58 kallisto_quant_19
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:57 kallisto_quant_19-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 16:59 kallisto_quant_19-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:09 kallisto_quant_2
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:02 kallisto_quant_20
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:02 kallisto_quant_20-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:04 kallisto_quant_20-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:05 kallisto_quant_21
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:07 kallisto_quant_28
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:07 kallisto_quant_29
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:15 kallisto_quant_3
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:09 kallisto_quant_30
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:11 kallisto_quant_31
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:11 kallisto_quant_37
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:13 kallisto_quant_38
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:13 kallisto_quant_39
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:21 kallisto_quant_4
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:15 kallisto_quant_40
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:17 kallisto_quant_41
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:17 kallisto_quant_47
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:19 kallisto_quant_48
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:19 kallisto_quant_49
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:32 kallisto_quant_5
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:22 kallisto_quant_50
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:24 kallisto_quant_57
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:24 kallisto_quant_57-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:25 kallisto_quant_57-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:27 kallisto_quant_58
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:27 kallisto_quant_58-G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:30 kallisto_quant_58-S
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:30 kallisto_quant_59
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:32 kallisto_quant_60
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:34 kallisto_quant_67
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:35 kallisto_quant_68
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:36 kallisto_quant_69
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:36 kallisto_quant_70
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:38 kallisto_quant_78
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:38 kallisto_quant_79
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:40 kallisto_quant_80
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:40 kallisto_quant_83
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:42 kallisto_quant_88
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:42 kallisto_quant_90
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:44 kallisto_quant_91
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:44 kallisto_quant_97
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:46 kallisto_quant_98
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:46 kallisto_quant_99
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:48 kallisto_quant_RESUB-116
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:48 kallisto_quant_RESUB-156
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:49 kallisto_quant_RESUB-36
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:50 kallisto_quant_RESUB-76
    drwxr-xr-x 2 shedurkin labmembers 4.0K Mar 18 17:51 kallisto_quant_RESUB-94
