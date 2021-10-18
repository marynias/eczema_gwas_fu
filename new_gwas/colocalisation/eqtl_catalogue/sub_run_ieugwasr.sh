#!/bin/bash
#SBATCH --job-name=ieugwasr
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


HOME=/mnt/storage/home/qh18484
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue


Rscript --vanilla $scripts/ieugwasr_coloc_local.R $my_interval $gene_list $gwas
