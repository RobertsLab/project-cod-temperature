#!/bin/bash

#SBATCH --job-name=subset-beagle
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/snp-testing-top-n-beagle-region.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 2-0:0:0

# Input files
base="/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2"
training="${base}/snp-training/regions"
testing="${base}/snp-testing/regions"
id_file="${base}/snp-testing/regions/test-IDs.txt"
output_beagle="testing.beagle.gz"

#### Now I need to filter my test-samples.beagle.gz and create 13 separate beagle files, one each containing the n sites/snps that
# I identified as predicting population:

numbers=(10 50 75 100 250 500 1000 1250 5000 10000 20000 35000 50000)

for n in "${numbers[@]}"
do

  # name of file to save
  sites=training.top_${n}_sites
  
  # filter sites 
  awk 'NR==FNR{c[$1]++;next};c[$1]' ${training}/${sites} <(zcat ${testing}/${output_beagle}) | gzip > ${testing}/testing-top-${n}.beagle.gz
 
done

## Separately create beagle file with all sites - if desired
#awk 'NR==FNR{c[$1]++;next};c[$1]' ${training}/training.all_sites <(zcat ${testing}/${output_beagle}) | gzip > ${testing}/testing-all-sites.beagle.gz

#for file in ${testing}/testing-top-*.gz; do echo "$file has $(zcat "$file" | wc -l) markers"; done

