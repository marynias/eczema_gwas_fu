#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/magma
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/magma
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Annotate the master file using MAGMA file generated for the pops pipeline.
Rscript --vanilla $scripts/magma_annotation.R $master/paternoster2015_master.csv \
$HOME/bin/pops/paternoster2015.genes.out \
paternoster2015_magma.csv
