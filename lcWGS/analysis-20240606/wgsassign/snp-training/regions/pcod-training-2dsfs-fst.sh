#!/bin/bash
#SBATCH -t 60-0:0:0
#SBATCH -p himem
#SBATCH --job-name=2dsfs
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/err-out/2dsfs-fst-regions_%A_%a.out
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --array=1-10

# Script to generate 2-dimentional site frequency spectrums, and Fst calculations among populations 
# Number in array = number of pairwise population combinations 

base=/home/lspencer/pcod-lcwgs-2023/analysis-20240606/wgsassign2/snp-training/regions
saf_dir=${base}/region-saf
sfs_dir=${base}/region-2dsfs
fst_dir=${base}/region-fst
angsd_dir=/opt/bioinformatics/bio/angsd/angsd-0.940/bin

# used this code to make pops-combos.txt file (all combinations of populations for pairwise 2dsfs and Fst); only need to run it once! 
#awk 'NR==FNR {a[NR]=$1; n=NR; next} {for (i=1; i<=n; i++) if (a[i] < $1) print a[i], $1}' OFS='\t' ${base}/regions.txt ${base}/regions.txt | sort -k1,1 > ${base}/regions-combos.txt

# Get the specific line for this array job
region1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $1}' ${base}/regions-combos.txt)
region2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $2}' ${base}/regions-combos.txt)

# Use the values of column1 and column2 as variables
echo "Calculating 2d-site frequency spectrum for $region1 x $region2"

outname=training.${region1}.${region2}

# get 2dsfs
${angsd_dir}/realSFS ${saf_dir}/wgsassign-training.${region1}.saf.idx ${saf_dir}/wgsassign-training.${region2}.saf.idx \
-maxIter 1000 -P 8 > ${sfs_dir}/${outname}

echo "Calculating Fst for $region1 x $region2"
 
# get fst
${angsd_dir}/realSFS fst index ${saf_dir}/wgsassign-training.${region1}.saf.idx ${saf_dir}/wgsassign-training.${region2}.saf.idx \
-sfs ${sfs_dir}/${outname} -whichFst 1 -fstout ${fst_dir}/${outname}

# Getting individual locus Fst stuff
${angsd_dir}/realSFS fst print ${fst_dir}/${outname}.fst.idx > ${fst_dir}/${outname}.fst

echo "Creating new Fst file containing only sites with Fst > 0"

# just report Fsts > 0 to minimize file size
awk '$4 != 0 {print $1, $2, $3/$4}' ${fst_dir}/${outname}.fst | awk '$3 > 0 {print}' | sort -k3 -rg > ${fst_dir}/${outname}.gt0.fst
