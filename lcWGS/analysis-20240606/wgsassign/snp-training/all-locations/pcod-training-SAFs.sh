#!/bin/bash
#SBATCH -t 60-0:0:0
#SBATCH -p himem
#SBATCH --job-name=wgsassign-saf
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/training-all-pops-saf.out
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL

module load bio/angsd/0.940

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606
training=${base}/wgsassign2/snp-training/all-locations

populations=$(cut -f5 ${training}/training_bams-list-pops.txt | sort | uniq)
sites=${base}/wgsassign2/pcod-refs_wholegenome_unlinked.sites
reference1=/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic.fa
reference2=/home/lspencer/references/pcod-ncbi/GCF_031168955.1_ASM3116895v1_genomic_rehead.fa #reheaded to have only chromosome ID, not other info
outdir=${training}/pop-saf

for pop in ${populations}
do
  grep -E ${pop} ${training}/training_bams-list-pops.txt | cut -f1 > ${training}/pop-saf/${pop}_bams-list.txt
  input_bams=${training}/pop-saf/${pop}_bams-list.txt
  outname=wgsassign-training.${pop}

  angsd -b ${input_bams} -out ${outdir}/${outname} \
  -gl 2 -domajorminor 3 -dosaf 1 -anc ${reference1} \
  -sites ${sites} -nThreads 8
done 
