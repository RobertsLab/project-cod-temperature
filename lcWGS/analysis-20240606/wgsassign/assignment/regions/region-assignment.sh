#!/bin/bash

#SBATCH --job-name=wgsassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/region-assign-top1250.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2
assign=${base}/assignment/regions
best_snps_beagle=${base}/snp-testing/regions/testing-top-1250.beagle.gz
best_snps_afs=${base}/snp-testing/regions/testing-LOOs/top-1250.pop_af.npy
exp_beagle_all=${base}/experimental.beagle.gz 
exp_beagle_best=${assign}/pcod-exp_top-1250.beagle
outname=${assign}/pcod-experimental-assign_top-1250-snps

# Create sites file containing the best SNPs to use for population assignment
zcat $best_snps_beagle | cut -f1 > ${assign}/reg-assign-1250.sites

# combine header row with filtered sites
zcat ${exp_beagle_all} | head -n 1 > ${exp_beagle_best}
awk 'NR==FNR{c[$1]++;next};c[$1]' ${assign}/reg-assign-1250.sites <(zcat ${exp_beagle_all}) >> ${exp_beagle_best}
gzip ${exp_beagle_best}

## IF BEAGLE FILE HAS NAMES OF SAMPLES IN HEADER ROW 
## Generate file listing samples in order they appear in beagle file using sample IDs in header row
#zcat ${exp_beagle_best}.gz | cut --complement -f1-3 | head -n 1 | tr '\t' '\n' | \
#sed -e 's/_AA//g' -e 's/_AB//g' -e 's/_BB//g' | uniq > ${assign}/exp-sample-order_reg.txt

# IF BEAGLE FILE ONLY HAS "Ind1", "Ind2" etc. in header row, copy the bams.list.txt file used in ANGSD to generated the experimental beagle file to get sample order 

# Run assignment to identify source population of experimental fish 
WGSassign --beagle ${exp_beagle_best}.gz --pop_af_file ${best_snps_afs} --get_pop_like --out ${outname} --threads 20
