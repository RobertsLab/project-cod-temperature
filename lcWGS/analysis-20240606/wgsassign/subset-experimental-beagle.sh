#!/bin/bash

#SBATCH --job-name=subset-beagle
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/subset-exp-from-merged-beagle.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -t 2-0:0:0

# Input files
base="/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2"
input_beagle="${base}/refs-exp-merged.beagle.gz"
cat ${base}/refs-exp-merged_samples-order.txt | grep "GM" > ${base}/exp_IDs.txt
id_file="${base}/exp_IDs.txt"
output_beagle="${base}/experimental.beagle.gz"

### First I need to subset my beagle file for only my test individuals 

# Extract IDs from the id_file into a regex pattern (exact match)
ids=$(awk '{print "^" $1 "$"}' $id_file | paste -sd '|' -)

# Process the input file
zcat $input_beagle | awk -v ids="$ids" '
BEGIN { FS=OFS="\t"; }
NR==1 {
    # Print the first three columns (marker, allele1, allele2)
    header_line = $1 OFS $2 OFS $3;
    for (i=4; i<=NF; i++) {
        # Remove suffixes and check if the column matches any ID
        col_name = gensub(/_AA$|_AB$|_BB$/, "", "g", $i)
        if (col_name ~ ids) {
            header_line = header_line OFS $i
            cols_to_print[i] = 1
        }
    }
    print header_line
}
NR>1 {
    line = $1 OFS $2 OFS $3;
    for (i=4; i<=NF; i++) {
        if (i in cols_to_print) {
            line = line OFS $i
        }
    }
    print line
}
' | gzip > $output_beagle
