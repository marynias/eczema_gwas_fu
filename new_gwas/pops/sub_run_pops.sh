#!/bin/bash

#SBATCH --job-name=pops
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=99:59:00
#SBATCH --mem=60G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

HOME=/mnt/storage/home/qh18484
pops=$HOME/scratch/pops

cd $pops

python pops.feature_selection.py \
	--features PoPS.features.txt.gz \
	--gene_results $gwas \
	--out $gwas 

for CHR in {1..22}
do
python pops.predict_scores.py \
	--gene_loc gene_loc.txt \
	--gene_results $gwas \
	--features PoPS.features.txt.gz \
	--selected_features ${gwas}.features\
	--control_features control.features \
	--chromosome $CHR \
	--out $gwas 
done