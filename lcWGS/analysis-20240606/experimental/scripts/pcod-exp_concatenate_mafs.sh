#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=pcod-exp_concat-mafs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/job_outfiles/pcod-exp_concatenate-mafs_%A.out

module unload bio/angsd/0.933
module load bio/angsd/0.933

for i in /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr1_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr2_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr3_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr4_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr5_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr6_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr7_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr8_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr9_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr10_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr11_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr12_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr13_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr14_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr15_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr16_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr17_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr18_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr19_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr20_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr21_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr22_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr23_wgassign.mafs.gz /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_Chr24_wgassign.mafs.gz
do zcat $i | tail -n +2 -q >> /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.mafs; done
cut -f 1,2,3,4 /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.mafs > /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.sites
gzip /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.mafs

angsd sites index /home/lspencer/pcod-lcwgs-2023/analysis-20240606/experimental/gls_wgassign/pcod-exp_wholegenome_wgassign.sites