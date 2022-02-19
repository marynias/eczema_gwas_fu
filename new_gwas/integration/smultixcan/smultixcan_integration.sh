#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/smultixcan
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/smultixcan
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
master=$HOME/scratch/new_gwas/integration/1_loci_prep
integration=$HOME/scratch/new_gwas/integration/final_integration

gwas_name="eczema21_discovery"

cd $analysis

#Annotate the master file.
Rscript --vanilla $scripts/smultixcan_annotation.R $master/${gwas_name}_master.csv \
shared_tissue_studies.csv \
${gwas_name}_smultixcan.csv

#Move final output file to the final integration folder.
cp ${gwas_name}_smultixcan.csv $integration