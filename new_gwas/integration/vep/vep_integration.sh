#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/vep
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/vep
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/scratch/new_gwas/genome
master=$HOME/scratch/new_gwas/integration/1_loci_prep
integration=$HOME/scratch/new_gwas/integration/final_integration

gwas_name="eczema21_discovery"

cd $analysis

#Annotate the master file with missense and intronic columns.
Rscript --vanilla $scripts/vepecho_annotation.R $master/${gwas_name}_master.csv \
ZIZHPgOugipEbG80.txt  \
${gwas_name}_vep.csv

#Move final output file to the final integration folder.
cp ${gwas_name}_vep.csv $integration
