#!/bin/bash

#SBATCH --job-name=snpsea-viz
#SBATCH --partition=mrcieu2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


snpsea=/mnt/storage/home/qh18484/scratch/snpsea

python $snpsea/bin/snpsea-barplot $output
#python $snpsea/bin/snpsea-heatmap $output
#Rscript $snpsea/bin/snpsea-type1error $output