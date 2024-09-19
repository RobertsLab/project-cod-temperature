#!/bin/bash
#SBATCH -t 60-0:0:0
#SBATCH -p himem
#SBATCH --job-name=wgsassign-saf
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/saf-region_%A_%a.out
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --array=1-5

module load bio/angsd/0.940

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2
training=${base}/snp-training/regions
mkdir ${training}/region-saf
outdir=${training}/region-saf
#cat training_bams-list_MR.txt | cut -f3 | sort | uniq > ${training}/regions.txt
regions_file=${training}/regions.txt
region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $regions_file)
sites=${base}/pcod-refs_wholegenome_unlinked.sites
reference1=/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa
reference2=/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_rehead.fa #reheaded to have only chromosome ID, not other info

grep -E $region ${training}/training_bams-list_MR.txt | cut -f1 > ${training}/region-saf/${region}_bams-list.txt
input_bams=${training}/region-saf/${region}_bams-list.txt
outname=wgsassign-training.${region}

angsd -b ${input_bams} -out ${outdir}/${outname} \
-gl 2 -domajorminor 3 -dosaf 1 -anc ${reference1} \
-sites ${sites} -nThreads 8
