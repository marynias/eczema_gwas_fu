#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/coloc
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/coloc
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Annotate the master file.
Rscript --vanilla $scripts/coloc_annotation.R $master/paternoster2015_master.csv \
gxp_coloc_all_results.txt \
paternoster2015_open_targets.csv
