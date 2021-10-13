#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/1_loci_prep
input=$HOME/scratch/new_gwas/loci_definition
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/1_loci_prep

gwas_name="eczema21_discovery"

cd $analysis

#Note: can use MendelVar to generate locus intervals, however if SNP is missing in the database,
#will have to find its proxy manually. May want to use the following synonym list:
$HOME/analysis/annotation/data_manipulation/rsid_synonyms.txt

#In that instance, MendelVar generated intervals for all the SNPs. The generated bed file is in:
cat $input/${gwas_name}_interval.bed

#First, make sure to run the scripts in $HOME/bin/eczema_gwas_fu/new_gwas/loci_definition.
#Create the basic table for each locus. Use only gene with assigned HGNC ID.
Rscript --vanilla $scripts/create_master_table.R $gwas/leadSNPs.${gwas_name}_cytoBand.bed \
$input/${gwas_name}_sorted_1Mbp_genes_processed ${gwas_name}_master.csv
