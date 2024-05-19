09-Hisat
================
Steven Roberts
15 May, 2024

- <a href="#1-reads" id="toc-1-reads">1 Reads</a>
- <a href="#2-genome" id="toc-2-genome">2 Genome</a>

# 1 Reads

``` bash
ls /home/shared/8TB_HDD_02/shedurkin/project-cod-temperature/output/05-cod-RNAseq-trimming/trimmed-reads
```

    100.flexbar_trim.R_1.fastq.gz
    100.flexbar_trim.R_2.fastq.gz
    100.flexbar_trim.R.log
    107.flexbar_trim.R_1.fastq.gz
    107.flexbar_trim.R_2.fastq.gz
    107.flexbar_trim.R.log
    108.flexbar_trim.R_1.fastq.gz
    108.flexbar_trim.R_2.fastq.gz
    108.flexbar_trim.R.log
    109.flexbar_trim.R_1.fastq.gz
    109.flexbar_trim.R_2.fastq.gz
    109.flexbar_trim.R.log
    10.flexbar_trim.R_1.fastq.gz
    10.flexbar_trim.R_2.fastq.gz
    10.flexbar_trim.R.log
    110.flexbar_trim.R_1.fastq.gz
    110.flexbar_trim.R_2.fastq.gz
    110.flexbar_trim.R.log
    117.flexbar_trim.R_1.fastq.gz
    117.flexbar_trim.R_2.fastq.gz
    117.flexbar_trim.R.log
    118.flexbar_trim.R_1.fastq.gz
    118.flexbar_trim.R_2.fastq.gz
    118.flexbar_trim.R.log
    119.flexbar_trim.R_1.fastq.gz
    119.flexbar_trim.R_2.fastq.gz
    119.flexbar_trim.R.log
    11.flexbar_trim.R_1.fastq.gz
    11.flexbar_trim.R_2.fastq.gz
    11.flexbar_trim.R.log
    120.flexbar_trim.R_1.fastq.gz
    120.flexbar_trim.R_2.fastq.gz
    120.flexbar_trim.R.log
    121.flexbar_trim.R_1.fastq.gz
    121.flexbar_trim.R_2.fastq.gz
    121.flexbar_trim.R.log
    127.flexbar_trim.R_1.fastq.gz
    127.flexbar_trim.R_2.fastq.gz
    127.flexbar_trim.R.log
    128.flexbar_trim.R_1.fastq.gz
    128.flexbar_trim.R_2.fastq.gz
    128.flexbar_trim.R.log
    129.flexbar_trim.R_1.fastq.gz
    129.flexbar_trim.R_2.fastq.gz
    129.flexbar_trim.R.log
    12.flexbar_trim.R_1.fastq.gz
    12.flexbar_trim.R_2.fastq.gz
    12.flexbar_trim.R.log
    131.flexbar_trim.R_1.fastq.gz
    131.flexbar_trim.R_2.fastq.gz
    131.flexbar_trim.R.log
    137.flexbar_trim.R_1.fastq.gz
    137.flexbar_trim.R_2.fastq.gz
    137.flexbar_trim.R.log
    138.flexbar_trim.R_1.fastq.gz
    138.flexbar_trim.R_2.fastq.gz
    138.flexbar_trim.R.log
    139.flexbar_trim.R_1.fastq.gz
    139.flexbar_trim.R_2.fastq.gz
    139.flexbar_trim.R.log
    13.flexbar_trim.R_1.fastq.gz
    13.flexbar_trim.R_2.fastq.gz
    13.flexbar_trim.R.log
    140.flexbar_trim.R_1.fastq.gz
    140.flexbar_trim.R_2.fastq.gz
    140.flexbar_trim.R.log
    147.flexbar_trim.R_1.fastq.gz
    147.flexbar_trim.R_2.fastq.gz
    147.flexbar_trim.R.log
    148.flexbar_trim.R_1.fastq.gz
    148.flexbar_trim.R_2.fastq.gz
    148.flexbar_trim.R.log
    149.flexbar_trim.R_1.fastq.gz
    149.flexbar_trim.R_2.fastq.gz
    149.flexbar_trim.R.log
    150.flexbar_trim.R_1.fastq.gz
    150.flexbar_trim.R_2.fastq.gz
    150.flexbar_trim.R.log
    18.flexbar_trim.R_1.fastq.gz
    18.flexbar_trim.R_2.fastq.gz
    18.flexbar_trim.R.log
    19.flexbar_trim.R_1.fastq.gz
    19.flexbar_trim.R_2.fastq.gz
    19.flexbar_trim.R.log
    19-G.flexbar_trim.R_1.fastq.gz
    19-G.flexbar_trim.R_2.fastq.gz
    19-G.flexbar_trim.R.log
    19-S.flexbar_trim.R_1.fastq.gz
    19-S.flexbar_trim.R_2.fastq.gz
    19-S.flexbar_trim.R.log
    1.flexbar_trim.R_1.fastq.gz
    1.flexbar_trim.R_2.fastq.gz
    1.flexbar_trim.R.log
    20.flexbar_trim.R_1.fastq.gz
    20.flexbar_trim.R_2.fastq.gz
    20.flexbar_trim.R.log
    20-G.flexbar_trim.R_1.fastq.gz
    20-G.flexbar_trim.R_2.fastq.gz
    20-G.flexbar_trim.R.log
    20-S.flexbar_trim.R_1.fastq.gz
    20-S.flexbar_trim.R_2.fastq.gz
    20-S.flexbar_trim.R.log
    21.flexbar_trim.R_1.fastq.gz
    21.flexbar_trim.R_2.fastq.gz
    21.flexbar_trim.R.log
    28.flexbar_trim.R_1.fastq.gz
    28.flexbar_trim.R_2.fastq.gz
    28.flexbar_trim.R.log
    29.flexbar_trim.R_1.fastq.gz
    29.flexbar_trim.R_2.fastq.gz
    29.flexbar_trim.R.log
    2.flexbar_trim.R_1.fastq.gz
    2.flexbar_trim.R_2.fastq.gz
    2.flexbar_trim.R.log
    30.flexbar_trim.R_1.fastq.gz
    30.flexbar_trim.R_2.fastq.gz
    30.flexbar_trim.R.log
    31.flexbar_trim.R_1.fastq.gz
    31.flexbar_trim.R_2.fastq.gz
    31.flexbar_trim.R.log
    37.flexbar_trim.R_1.fastq.gz
    37.flexbar_trim.R_2.fastq.gz
    37.flexbar_trim.R.log
    38.flexbar_trim.R_1.fastq.gz
    38.flexbar_trim.R_2.fastq.gz
    38.flexbar_trim.R.log
    39.flexbar_trim.R_1.fastq.gz
    39.flexbar_trim.R_2.fastq.gz
    39.flexbar_trim.R.log
    3.flexbar_trim.R_1.fastq.gz
    3.flexbar_trim.R_2.fastq.gz
    3.flexbar_trim.R.log
    40.flexbar_trim.R_1.fastq.gz
    40.flexbar_trim.R_2.fastq.gz
    40.flexbar_trim.R.log
    41.flexbar_trim.R_1.fastq.gz
    41.flexbar_trim.R_2.fastq.gz
    41.flexbar_trim.R.log
    47.flexbar_trim.R_1.fastq.gz
    47.flexbar_trim.R_2.fastq.gz
    47.flexbar_trim.R.log
    48.flexbar_trim.R_1.fastq.gz
    48.flexbar_trim.R_2.fastq.gz
    48.flexbar_trim.R.log
    49.flexbar_trim.R_1.fastq.gz
    49.flexbar_trim.R_2.fastq.gz
    49.flexbar_trim.R.log
    4.flexbar_trim.R_1.fastq.gz
    4.flexbar_trim.R_2.fastq.gz
    4.flexbar_trim.R.log
    50.flexbar_trim.R_1.fastq.gz
    50.flexbar_trim.R_2.fastq.gz
    50.flexbar_trim.R.log
    57.flexbar_trim.R_1.fastq.gz
    57.flexbar_trim.R_2.fastq.gz
    57.flexbar_trim.R.log
    57-G.flexbar_trim.R_1.fastq.gz
    57-G.flexbar_trim.R_2.fastq.gz
    57-G.flexbar_trim.R.log
    57-S.flexbar_trim.R_1.fastq.gz
    57-S.flexbar_trim.R_2.fastq.gz
    57-S.flexbar_trim.R.log
    58.flexbar_trim.R_1.fastq.gz
    58.flexbar_trim.R_2.fastq.gz
    58.flexbar_trim.R.log
    58-G.flexbar_trim.R_1.fastq.gz
    58-G.flexbar_trim.R_2.fastq.gz
    58-G.flexbar_trim.R.log
    58-S.flexbar_trim.R_1.fastq.gz
    58-S.flexbar_trim.R_2.fastq.gz
    58-S.flexbar_trim.R.log
    59.flexbar_trim.R_1.fastq.gz
    59.flexbar_trim.R_2.fastq.gz
    59.flexbar_trim.R.log
    5.flexbar_trim.R_1.fastq.gz
    5.flexbar_trim.R_2.fastq.gz
    5.flexbar_trim.R.log
    60.flexbar_trim.R_1.fastq.gz
    60.flexbar_trim.R_2.fastq.gz
    60.flexbar_trim.R.log
    67.flexbar_trim.R_1.fastq.gz
    67.flexbar_trim.R_2.fastq.gz
    67.flexbar_trim.R.log
    68.flexbar_trim.R_1.fastq.gz
    68.flexbar_trim.R_2.fastq.gz
    68.flexbar_trim.R.log
    69.flexbar_trim.R_1.fastq.gz
    69.flexbar_trim.R_2.fastq.gz
    69.flexbar_trim.R.log
    70.flexbar_trim.R_1.fastq.gz
    70.flexbar_trim.R_2.fastq.gz
    70.flexbar_trim.R.log
    78.flexbar_trim.R_1.fastq.gz
    78.flexbar_trim.R_2.fastq.gz
    78.flexbar_trim.R.log
    79.flexbar_trim.R_1.fastq.gz
    79.flexbar_trim.R_2.fastq.gz
    79.flexbar_trim.R.log
    80.flexbar_trim.R_1.fastq.gz
    80.flexbar_trim.R_2.fastq.gz
    80.flexbar_trim.R.log
    83.flexbar_trim.R_1.fastq.gz
    83.flexbar_trim.R_2.fastq.gz
    83.flexbar_trim.R.log
    88.flexbar_trim.R_1.fastq.gz
    88.flexbar_trim.R_2.fastq.gz
    88.flexbar_trim.R.log
    90.flexbar_trim.R_1.fastq.gz
    90.flexbar_trim.R_2.fastq.gz
    90.flexbar_trim.R.log
    91.flexbar_trim.R_1.fastq.gz
    91.flexbar_trim.R_2.fastq.gz
    91.flexbar_trim.R.log
    97.flexbar_trim.R_1.fastq.gz
    97.flexbar_trim.R_2.fastq.gz
    97.flexbar_trim.R.log
    98.flexbar_trim.R_1.fastq.gz
    98.flexbar_trim.R_2.fastq.gz
    98.flexbar_trim.R.log
    99.flexbar_trim.R_1.fastq.gz
    99.flexbar_trim.R_2.fastq.gz
    99.flexbar_trim.R.log
    RESUB-116.flexbar_trim.R_1.fastq.gz
    RESUB-116.flexbar_trim.R_2.fastq.gz
    RESUB-116.flexbar_trim.R.log
    RESUB-156.flexbar_trim.R_1.fastq.gz
    RESUB-156.flexbar_trim.R_2.fastq.gz
    RESUB-156.flexbar_trim.R.log
    RESUB-36.flexbar_trim.R_1.fastq.gz
    RESUB-36.flexbar_trim.R_2.fastq.gz
    RESUB-36.flexbar_trim.R.log
    RESUB-76.flexbar_trim.R_1.fastq.gz
    RESUB-76.flexbar_trim.R_2.fastq.gz
    RESUB-76.flexbar_trim.R.log
    RESUB-94.flexbar_trim.R_1.fastq.gz
    RESUB-94.flexbar_trim.R_2.fastq.gz
    RESUB-94.flexbar_trim.R.log
    trimmed_fastq_checksums.md5

Need to find treatment file…

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
