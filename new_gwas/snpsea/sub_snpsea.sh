#!/bin/bash

#SBATCH --job-name=snpsea
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

snpsea=/mnt/storage/home/qh18484/scratch/snpsea

$snpsea/bin/snpsea --snps $my_gwas \
--gene-matrix  $gene_matrix \
--gene-intervals  $snpsea/NCBIgenes2013.bed.gz \
--snp-intervals  $snpsea/TGP2011.bed.gz \
--null-snps $snpsea/Lango2010.txt.gz \
--out $output \
--slop 10e3 \
--threads 2 \
--null-snpsets 100 \
--min-observations 100 \
--max-iterations 1e7 