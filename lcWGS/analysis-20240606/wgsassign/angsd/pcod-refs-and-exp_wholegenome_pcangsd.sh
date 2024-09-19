#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=pcangsd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign/err-out/pcangsd_all_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/refs-exp-merged.beagle.gz -o /home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/angsd/refs-exp-merged_pca -sites_save -pcadapt

# Save order of samples to add column and row names 
zcat refs-exp-merged.beagle.gz | head -n 1 |  tr '\t' '\n' | tail -n +4 | sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > refs-exp-merged_samples-order.txt
