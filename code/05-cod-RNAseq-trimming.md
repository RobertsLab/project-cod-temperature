05-cod-RNAseq-trimming
================
Kathleen Durkin
2024-03-04

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

- Trimmed reads: `*flexbar_trim.25bp.fastq.gz`

(info about library prep/sequencing)

------------------------------------------------------------------------

# Create a Bash variables file

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

# Download raw sRNAseq reads

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

## Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

cd "${raw_reads_dir}"

md5sum checksums.md5 --check
```

# FastQC/MultiQC on raw reads

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

# View directory contents
ls -lh ${raw_fastqc_dir}
```

######## 

Not sure yet if I’m going to trim, ignore for now \########

# Create adapters FastA for use with [flexbar](https://github.com/seqan/flexbar) trimming

``` bash
# Load bash variables into memory
source .bashvars

echo "Creating adapters FastA."
echo ""
adapter_count=0

# Check for adapters file first
# Then create adapters file if doesn't exist
if [ -f "${output_dir_top}/${NEB_adapters_fasta}" ]; then
  echo "${output_dir_top}/${NEB_adapters_fasta} already exists. Nothing to do."
else
  for adapter in "${first_adapter}" "${second_adapter}"
  do
    adapter_count=$((adapter_count + 1))
    printf ">%s\n%s\n" "adapter_${adapter_count}" "${adapter}"
  done >> "${output_dir_top}/${NEB_adapters_fasta}"
fi

echo ""
echo "Adapters FastA:"
echo ""
cat "${output_dir_top}/${NEB_adapters_fasta}"
echo ""
```

# Concatenate reads (if run on multiple lanes)

# Trimming with [flexbar](https://github.com/seqan/flexbar)

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
# Uses NEB adapter file
# --adapter-pair-overlap ON: Recommended by NEB sRNA kit
# --qtrim-threshold 25: Minimum quality
# --qtrim-format i1.8: Sets sequencer as illumina
# --post-trim-length: Trim reads from 3' end to max length
# --target: Sets file naming patterns
# --zip-output GZ: Sets type of compression. GZ = gzip

# --adapter-min-overlap 7: Minimum overlap between adapter and read for trimming consideration
# --adapter-trim-end RIGHT: Trim adapters from the right (3') end
# --adapter-trim-position TAIL: Adapters should be found at the end (tail) of the reads
# --adapter-min-length 7: Sets the minimum length of an adapter. Adapters shorter than this length will be ignored
# --adapter-seq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA: Specifies the first Illumina adapter sequence to be trimmed
# --adapter-seq AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT: Specifies the second Illumina adapter sequence to be trimmed
# --qtrim-format sanger: Specifies the quality score format. In this case, it's set to Sanger format (Phred+33)
# --qtrim-threshold 20: Sets the quality threshold for trimming (bases w/ score < 20 will be trimmed)
# --qtrim-position BOTH: Quality trimming should be performed on both ends of the reads
# --qtrim-final: Quality trimming should be performed as the final step in the trimming process
# --min-read-length 50: The minimum length a read must have after trimming; reads shorter than this length will be discarded
# --threads 20: Number of threads (CPU cores) to use for parallel processing
# --pre-trim-left 0 and --pre-trim-right 0: Specifies the number of bases to be removed from the 5' (left) and 3' (right) ends of the reads before any adapter or quality trimming. In this case, no bases are removed.
# --pre-trim-phred 33: Specifies the Phred quality score offset used in the input files. For Illumina data, it's typically 33.
#--input1 and --input2: Input FASTQ files for the paired-end reads.
#--output1 and --output2: Output files for the trimmed paired-end reads.

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
  --adapters ${output_dir_top}/${NEB_adapters_fasta} \
  --adapter-pair-overlap ON \
  --qtrim-format i1.8 \
  --qtrim-threshold 25 \
  --post-trim-length ${max_read_length} \
  --threads ${threads} \
  --target "${trimmed_reads_dir}/${R1_sample_name%%_*}.flexbar_trim.${max_read_length}bp" \
  --zip-output GZ
        
    # Move to trimmed directory
    # This is done so checksums file doesn't include excess path
    cd ${trimmed_reads_dir}

    # Generate md5 checksums for newly trimmed files
    {
      md5sum "${R1_sample_name%%_*}.flexbar_trim.${max_read_length}bp_1.fastq.gz"
      md5sum "${R2_sample_name%%_*}.flexbar_trim.${max_read_length}bp_2.fastq.gz"
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

# FastQC/MultiQC on trimmed reads

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

# View directory contents
ls -lh ${trimmed_fastqc_dir}
```

# Summary