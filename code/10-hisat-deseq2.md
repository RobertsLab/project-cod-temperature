10-Hisat
================
Steven Roberts
19 May, 2024

- <a href="#1-reads" id="toc-1-reads">1 Reads</a>
  - <a href="#11-multiqc" id="toc-11-multiqc">1.1 multiqc</a>
- <a href="#2-genome" id="toc-2-genome">2 Genome</a>
- <a href="#3-hisat" id="toc-3-hisat">3 Hisat</a>

Alt splice test run

# 1 Reads

``` bash
ls ../data/reads/*
```

    ../data/reads/100.trimmed.R1.fastq.gz
    ../data/reads/100.trimmed.R2.fastq.gz
    ../data/reads/107.trimmed.R1.fastq.gz
    ../data/reads/107.trimmed.R2.fastq.gz
    ../data/reads/108.trimmed.R1.fastq.gz
    ../data/reads/108.trimmed.R2.fastq.gz
    ../data/reads/109.trimmed.R1.fastq.gz
    ../data/reads/109.trimmed.R2.fastq.gz
    ../data/reads/10.trimmed.R1.fastq.gz
    ../data/reads/10.trimmed.R2.fastq.gz
    ../data/reads/110.trimmed.R1.fastq.gz
    ../data/reads/110.trimmed.R2.fastq.gz
    ../data/reads/116.trimmed.R1.fastq.gz
    ../data/reads/116.trimmed.R2.fastq.gz
    ../data/reads/11.trimmed.R1.fastq.gz
    ../data/reads/11.trimmed.R2.fastq.gz
    ../data/reads/12.trimmed.R1.fastq.gz
    ../data/reads/12.trimmed.R2.fastq.gz
    ../data/reads/13.trimmed.R1.fastq.gz
    ../data/reads/13.trimmed.R2.fastq.gz
    ../data/reads/18.trimmed.R1.fastq.gz
    ../data/reads/18.trimmed.R2.fastq.gz
    ../data/reads/19.trimmed.R1.fastq.gz
    ../data/reads/19.trimmed.R2.fastq.gz
    ../data/reads/1.trimmed.R1.fastq.gz
    ../data/reads/1.trimmed.R2.fastq.gz
    ../data/reads/20.trimmed.R1.fastq.gz
    ../data/reads/20.trimmed.R2.fastq.gz
    ../data/reads/21.trimmed.R1.fastq.gz
    ../data/reads/21.trimmed.R2.fastq.gz
    ../data/reads/28.trimmed.R1.fastq.gz
    ../data/reads/28.trimmed.R2.fastq.gz
    ../data/reads/29.trimmed.R1.fastq.gz
    ../data/reads/29.trimmed.R2.fastq.gz
    ../data/reads/2.trimmed.R1.fastq.gz
    ../data/reads/2.trimmed.R2.fastq.gz
    ../data/reads/30.trimmed.R1.fastq.gz
    ../data/reads/30.trimmed.R2.fastq.gz
    ../data/reads/31.trimmed.R1.fastq.gz
    ../data/reads/31.trimmed.R2.fastq.gz
    ../data/reads/36.trimmed.R1.fastq.gz
    ../data/reads/36.trimmed.R2.fastq.gz
    ../data/reads/3.trimmed.R1.fastq.gz
    ../data/reads/3.trimmed.R2.fastq.gz
    ../data/reads/4.trimmed.R1.fastq.gz
    ../data/reads/4.trimmed.R2.fastq.gz
    ../data/reads/5.trimmed.R1.fastq.gz
    ../data/reads/5.trimmed.R2.fastq.gz
    ../data/reads/78.trimmed.R1.fastq.gz
    ../data/reads/78.trimmed.R2.fastq.gz
    ../data/reads/79.trimmed.R1.fastq.gz
    ../data/reads/79.trimmed.R2.fastq.gz
    ../data/reads/80.trimmed.R1.fastq.gz
    ../data/reads/80.trimmed.R2.fastq.gz
    ../data/reads/83.trimmed.R1.fastq.gz
    ../data/reads/83.trimmed.R2.fastq.gz
    ../data/reads/88.trimmed.R1.fastq.gz
    ../data/reads/88.trimmed.R2.fastq.gz
    ../data/reads/90.trimmed.R1.fastq.gz
    ../data/reads/90.trimmed.R2.fastq.gz
    ../data/reads/91.trimmed.R1.fastq.gz
    ../data/reads/91.trimmed.R2.fastq.gz
    ../data/reads/92.trimmed.R1.fastq.gz
    ../data/reads/92.trimmed.R2.fastq.gz
    ../data/reads/94.trimmed.R1.fastq.gz
    ../data/reads/94.trimmed.R2.fastq.gz
    ../data/reads/97.trimmed.R1.fastq.gz
    ../data/reads/97.trimmed.R2.fastq.gz
    ../data/reads/98.trimmed.R1.fastq.gz
    ../data/reads/98.trimmed.R2.fastq.gz
    ../data/reads/99.trimmed.R1.fastq.gz
    ../data/reads/99.trimmed.R2.fastq.gz
    ../data/reads/splice-test-files.txt

Need to find treatment file…

## 1.1 multiqc

``` bash
/home/sam/programs/mambaforge/bin/multiqc \
../data/reads/*gz
```

# 2 Genome

<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_031168955.1/>

![genome](http://gannet.fish.washington.edu/seashell/snaps/Monosnap_Gadus_macrocephalus_genome_assembly_ASM3116895v1_-_NCBI_-_NLM_2024-05-15_09-20-59.png)

``` bash
cd ../data

/home/shared/datasets download genome accession GCF_031168955.1 --include gff3,gtf,rna,cds,protein,genome,seq-report
```

``` bash
cd ../data 
unzip ncbi_dataset.zip
```

``` bash
ls /home/shared/8TB_HDD_03/sr320/github/project-cod-temperature/data/ncbi_dataset/data/GCF_031168955.1
```

    cds_from_genomic.fna
    GCF_031168955.1_ASM3116895v1_genomic.fna
    genomic.gff
    genomic.gtf
    protein.faa
    rna.fna
    sequence_report.jsonl

# 3 Hisat

``` bash
/home/shared/hisat2-2.2.1/hisat2_extract_exons.py \
../data/ncbi_dataset/data/GCF_031168955.1/genomic.gtf \
> ../output/10-hisat-deseq2/m_exon.tab
```

``` bash
/home/shared/hisat2-2.2.1/hisat2_extract_splice_sites.py \
../data/ncbi_dataset/data/GCF_031168955.1/genomic.gtf \
> ../output/10-hisat-deseq2/m_spice_sites.tab
```

``` bash
echo "10-hisat-deseq2/GCF*" >> ../output/.gitignore
```

``` bash
/home/shared/hisat2-2.2.1/hisat2-build \
../data/ncbi_dataset/data/GCF_031168955.1/GCF_031168955.1_ASM3116895v1_genomic.fna \
../output/10-hisat-deseq2/GCF_031168955.1.index \
--exon ../output/10-hisat-deseq2/m_exon.tab \
--ss ../output/10-hisat-deseq2/m_spice_sites.tab \
-p 20 \
../data/ncbi_dataset/data/GCF_031168955.1/genomic.gtf \
2> ../output/10-hisat-deseq2/hisat2-build_stats.txt
```

``` bash
echo "10-hisat-deseq2/*sam" >> ../output/.gitignore
```

``` bash
find ../data/reads/*.trimmed.R1.fastq.gz \
| xargs basename -s .trimmed.R1.fastq.gz | xargs -I{} \
/home/shared/hisat2-2.2.1/hisat2 \
-x ../output/10-hisat-deseq2/GCF_031168955.1.index \
--dta \
-p 20 \
-1 ../data/reads/{}.trimmed.R1.fastq.gz \
-2 ../data/reads/{}.trimmed.R2.fastq.gz \
-S ../output/10-hisat-deseq2/{}.sam \
2> ../output/10-hisat-deseq2/hisat.out
```

``` bash
for samfile in ../output/10-hisat-deseq2/*.sam; do
  bamfile="${samfile%.sam}.bam"
  sorted_bamfile="${samfile%.sam}.sorted.bam"
  /home/shared/samtools-1.12/samtools view -bS -@ 20 "$samfile" > "$bamfile"
  /home/shared/samtools-1.12/samtools sort -@ 20 "$bamfile" -o "$sorted_bamfile"
  /home/shared/samtools-1.12/samtools index -@ 20 "$sorted_bamfile"
done
```
