#!/bin/bash

#SBATCH --job-name=snp-testing
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/snp-testing-loos-regions.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -t 2-0:0:0

source ~/.bashrc
mamba activate WGSassign

mkdir testing-LOOs
base1=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/snp-testing
base2=$base1/regions
ids=${base2}/test-IDs.txt
numbers=(10 50 75 100 250 500 1000 1250 5000 10000 20000 35000 50000)

for n in "${numbers[@]}"
do

	# ex. amre.testing.ind85.ds_2x.sites-filter.top_100000_each.beagle.gz
	input_beagle=${base2}/testing-top-${n}.beagle.gz
	outname=${base2}/testing-LOOs/top-${n}

	# Get likelihoods for leave-one-out assignment within known reference populations
	# Output = 1) reference.popAF.npy, 2) reference.pop_like_LOO.txt
	WGSassign --beagle ${input_beagle} --pop_af_IDs ${ids} \
	--get_reference_af --loo --out ${outname} --threads 20

done

## Separately run this on all sites (not part of for loop), if desired
#WGSassign --beagle ${base1}/testing.beagle.gz --pop_af_IDs ${ids} \
#--get_reference_af --loo --out ${base2}/testing-LOOs/all-sites --threads 20
