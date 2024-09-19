#!/bin/bash

#SBATCH --job-name=wgsassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/pop-assign-top50.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2
best_snps_beagle=${base}/snp-testing/all-locations/testing-top-50.beagle.gz
sites=${base}/assignment/all-locations/pop-assign-50.sites
best_snps_afs=${base}/snp-testing/all-locations/testing-LOOs/top-50.pop_af.npy
exp_beagle_all=${base}/experimental.beagle.gz
exp_beagle_best=${base}/assignment/all-locations/pcod-exp_top-50.beagle
outname=${base}/assignment/all-locations/pcod-experimental-assign_top-50-snps

# Create sites file containing the best SNPs to use for population assignment
zcat $best_snps_beagle | cut -f1 > ${sites}

# combine header row with filtered sites
zcat ${exp_beagle_all} | head -n 1 > ${exp_beagle_best}
awk 'NR==FNR{c[$1]++;next};c[$1]' ${sites} <(zcat ${exp_beagle_all}) >> ${exp_beagle_best}
gzip ${exp_beagle_best}

## IF BEAGLE FILE HAS NAMES OF SAMPLES IN HEADER ROW 
## Generate file listing samples in order they appear in beagle file using sample IDs in header row
#zcat ${exp_beagle_best}.gz | cut --complement -f1-3 | head -n 1 | tr '\t' '\n' | \
#sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > ${base}/assignment/exp-sample-order_pop.txt

## IF BEAGLE FILE ONLY HAS "Ind1", "Ind2" etc. in header row, copy the bams.list.txt file used in ANGSD to generated the experimental beagle file to get sample order 

# Run assignment to identify source population of experimental fish 
WGSassign --beagle ${exp_beagle_best}.gz --pop_af_file ${best_snps_afs} --get_pop_like --out ${outname} --threads 20
