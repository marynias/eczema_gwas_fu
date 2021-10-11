#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/mendelvar
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/mendelvar
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Need to add a column for genes within MendelVar_enriched_pathway/phenotypes and MendelVar_skin_associated phenotypes (DO, HP keywords: skin, derma, kera). Run LD-based interval of r2 > 0.8 or 0.6 depending on which better enrichment results.* Remember to remove asterisks from output for sorting. For skin-related phenotypes consider all genes within 1 Mbp interval.

#Also have a glance at disease_overlap file and look for any additional genes with descriptions with skin, derma, kera keywords.

#*if lead SNP is missing in the database, will have to find its proxy manually and generate an interval that way.
#May want to use the following synonym list:
ls $HOME/analysis/annotation/data_manipulation/rsid_synonyms.txt

#Start off by coping the "parsed" INRICH output from MendelVar into the directory. Make sure to remove all the asterisk in the results.
Rscript --vanilla $scripts/mendelvar_annotation.R $master/paternoster2015_master.csv \
 paternoster2015_mendelvar_significant.csv \
 paternoster2015_mendelvar.csv

 