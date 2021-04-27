#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=ucgd-owner-guest
#SBATCH --partition=ucgd-owner-guest
#SBATCH -o get_variants.%j.log

## Define variables
dir=/uufs/chpc.utah.edu/common/home/u1240855/git/somccr

## Load modules
module use $HOME/MyModules
module load miniconda3/latest

## Remove existing files
echo "Removing preexisting bed files"
rm $dir/data/output/somccr_v3/tcga.*.bed
rm $dir/data/output/somccr_v3/icgc.*.bed

## Process TCGA variants
echo "Processing TCGA WES variants"
python $dir/script/somccr_v3/python/get_tcga.py $dir/data/output/somccr_v3/tcga.$SLURM_JOB_ID.bed.gz

## Process ICGC variants
echo "Processing ICGC WGS variants"
python $dir/script/somccr_v3/python/get_icgc.py $dir/data/output/somccr_v3/icgc.$SLURM_JOB_ID.bed.gz

