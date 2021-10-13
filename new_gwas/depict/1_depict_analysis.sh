#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/depict/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/depict
gwas=$HOME/scratch/new_gwas/gwas/raw

gwas_name="eczema21_discovery"

cd $analysis

#####Top independent loci analysis
#Top p-value 
cat $gwas/leadSNPs.${gwas_name}.index | cut -f3 >${gwas_name}_depict.top

#Before every run change parameters within the depict.py file in /mnt/storage/home/qh18484/scratch/depict
sbatch $scripts/sub_run_depict.sh