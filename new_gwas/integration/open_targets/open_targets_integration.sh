#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/open_targets
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/open_targets
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
master=$HOME/scratch/new_gwas/integration/1_loci_prep
results=$HOME/scratch/new_gwas/open_targets
integration=$HOME/scratch/new_gwas/integration/final_integration

gwas_name="eczema21_discovery"

cd $analysis

#Annotate the master file.
Rscript --vanilla $scripts/open_targets_annotation.R $master/${gwas_name}_master.csv \
$results/${gwas_name}_open_targets.txt \
${gwas_name}_open_targets.csv

#Move final output file to the final integration folder.
cp ${gwas_name}_open_targets.csv $integration
