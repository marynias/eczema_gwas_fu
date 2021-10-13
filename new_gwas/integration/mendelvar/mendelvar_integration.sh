#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/mendelvar
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/mendelvar
master=$HOME/scratch/new_gwas/integration/1_loci_prep

gwas_name="eczema21_discovery"

cd $analysis

#Need to add a column for genes within MendelVar_enriched_pathway/phenotypes 
#and MendelVar_skin_associated phenotypes (DO, HP keywords: skin, derma, kera). 
#Remember to remove asterisks from output for sorting. For skin-related phenotypes consider all genes within 1 Mbp interval.

for a in *parsed
do 
sed -i 's/[*]\+//g' $a
done


#*if lead SNP is missing in the database, will have to find its proxy manually and generate an interval that way.
#May want to use the following synonym list: Not the case here
ls $HOME/analysis/annotation/data_manipulation/rsid_synonyms.txt

#Start off by coping the "parsed" INRICH output from MendelVar into the directory. Make sure to remove all the asterisk in the results.
Rscript --vanilla $scripts/mendelvar_annotation.R $master/${gwas_name}_master.csv \
 ${gwas_name}_mendelvar_pathways_significant.csv \
 ${gwas_name}_mendelvar.csv

 