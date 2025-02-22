#!/bin/bash
#SBATCH -t 60-0:0:0
#SBATCH --job-name=top-snps
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/pull-top-snps-pops.out
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/snp-training/all-locations
numbers=(10 50 100 500 1000 5000 10000 15000 20000 25000 35000 45000)

for n in "${numbers[@]}"
do

  # name of file to save 
  name=training.top_${n}_sites

  # Print the top n most differentiated sites from all pairwise population comparisons, pull chrom+site, sort and pull unique ones. 
  head -n ${n} ${base}/pop-fst/*.gt0.fst | cut -d " " -f1-2 | grep "NC" | sort | uniq > ${base}/${name}.tmp

  # Replace space between chromasome and location with a "_" (to match beagle files) & remove temp. file
  tr -s ' ' '_' < ${base}/${name}.tmp > ${base}/${name}
  rm ${base}/${name}.tmp

done

cat /home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/pcod-refs_wholegenome_unlinked.sites | \
tr -s ' ' '_' <  > ${base}/training.all_sites
