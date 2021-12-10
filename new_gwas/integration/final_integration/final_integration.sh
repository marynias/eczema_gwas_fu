#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/final_integration
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/final_integration
master=$HOME/scratch/new_gwas/integration/1_loci_prep
coloc=$HOME/scratch/new_gwas/colocalisation/eqtl_catalogue
gwas=$HOME/scratch/new_gwas/gwas/raw

cd $analysis

gwas_name="eczema21_discovery"
#Integrate all individual results tables. At the end may want to split the table by rsid or locus.
#Check that all the tables have the same number of rows and then join by rsid,gene name,cytoband etc.

Rscript --vanilla $scripts/create_summary_tables.R $master/${gwas_name}_master.csv \
$gwas/Table_signHits_7th_500kb_run_updated_RSID.csv \
$gwas_name 


##Plot the overall results
Rscript --vanilla $scripts/plot_summary_table_ranking.R $gwas_name

##Plot colocalisation results.
Rscript --vanilla $scripts/plot_summary_table_coloc.R $coloc/${gwas_name}_all_coloc.txt \
$gwas/leadSNPs.${gwas_name}_cytoBand.bed $gwas_name

#Plots using manually chosen top 1 gene at a locus in case of ties.
Rscript --vanilla $scripts/create_summary_tables_chosenhits.R $master/${gwas_name}_master.csv \
$gwas/Table_signHits_7th_500kb_run_updated_RSID.csv \
$gwas_name ${gwas_name}_top_genes.txt 

##Plot the overall results
Rscript --vanilla $scripts/plot_summary_table_ranking_chosenhits.R $gwas_name

