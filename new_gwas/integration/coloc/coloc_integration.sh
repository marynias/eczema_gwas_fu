#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/coloc
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/coloc
master=$HOME/scratch/new_gwas/integration/1_loci_prep
coloc=$HOME/scratch/new_gwas/colocalisation/eqtl_catalogue
integration=$HOME/scratch/new_gwas/integration/final_integration

gwas_name="eczema21_discovery"

cd $analysis

#Annotate the master file.
Rscript --vanilla $scripts/coloc_annotation.R $master/${gwas_name}_master.csv \
$coloc/${gwas_name}_all_coloc.txt \
${gwas_name}_coloc.csv

#Move final output file to the final integration folder.
cp ${gwas_name}_coloc.csv $integration
