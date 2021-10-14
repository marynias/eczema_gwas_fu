#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/final_integration
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/final_integration

cd $analysis

#Integrate all individual results tables. At the end may want to split the table by rsid or locus.
#Check that all the tables have the same number of rows and then join by rsid,gene name,cytoband etc.

Rscript --vanilla $scripts/create_summary_tables.R

##Plot the overall results
Rsript --vanilla $scripts/plot_summary_table_ranking.R

##Plot colocalisation results.
Rsript --vanilla $scripts/plot_summary_table_coloc.R