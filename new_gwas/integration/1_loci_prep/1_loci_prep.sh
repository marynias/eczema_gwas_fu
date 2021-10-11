#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/1_loci_prep
input=$HOME/new_gwas/loci_definition/eqtl_catalogue
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/1_loci_prep


cd $analysis

#Note: can use MendelVar to generate locus intervals, however if SNP is missing in the database,
#will have to find its proxy manually. May want to use the following synonym list:
$HOME/analysis/annotation/data_manipulation/rsid_synonyms.txt

#First, make sure to run the script in $HOME/bin/eczema_gwas_fu/new_gwas/loci_definition.

#Create the basic table for each locus. Use only gene with assigned HGNC ID.
Rscript --vanilla $scripts/create_master_table.R $input/paternoster_2015_index_snps_cytoband.bed \
$input/paternoster_2015_index_snps_sorted_1Mbp_genes_processed.bed paternoster2015_master.csv
