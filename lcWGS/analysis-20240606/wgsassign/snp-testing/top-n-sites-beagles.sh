#!/bin/bash

#SBATCH --job-name=subset-beagle
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/snp-testing-subset-beagle.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 2-0:0:0

# Input files
base="/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2"
id_file="${base}/snp-testing/test-IDs.txt"
output_beagle="${base}/snp-testing/testing.beagle.gz"

#### Now I need to filter my test-samples.beagle.gz and create 12 separate beagle files, one each containing the n sites/snps that
# I identified as predicting population:

numbers=(10 50 100 500 1000 5000 10000 15000 20000 25000 35000 45000)

for n in "${numbers[@]}"
do

  # name of file to save
  sites=training.top_${n}_sites
  
  # filter sites 
  awk 'NR==FNR{c[$1]++;next};c[$1]' ${base}/snp-training/all-locations/${sites} <(zcat ${output_beagle}) | gzip > ${base}/snp-testing/all-locations/testing-top-${n}.beagle.gz
 
done

# Separately create beagle file with all sites
awk 'NR==FNR{c[$1]++;next};c[$1]' ${base}/snp-training/all-locations/training.all_sites <(zcat ${input_beagle}) | gzip > ${base}/snp-testing/all-locations/testing-all-sites.beagle.gz

for file in ${base}/snp-testing/all-locations/testing-top-*.gz; do echo "$file has $(zcat "$file" | wc -l) markers"; done

