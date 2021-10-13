#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/magma
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/magma
master=$HOME/scratch/new_gwas/integration/1_loci_prep

gwas_name="eczema21_discovery"

cd $analysis

#Annotate the master file using MAGMA file generated for the pops pipeline.
Rscript --vanilla $scripts/magma_annotation.R $master/${gwas_name}_master.csv \
$HOME/scratch/pops/${gwas_name}.genes.out \
${gwas_name}_magma.csv
