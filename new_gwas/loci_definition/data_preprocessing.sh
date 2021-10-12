#!/bin/bash
HOME=/mnt/storage/home/qh18484
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/loci_definition
gwas=$HOME/scratch/new_gwas/gwas/raw

cd $gwas

gwas_name="eczema21_discovery"

#Read the summary stats into R, rename columns so that consistent with Paternoster2015 files and output as TSV.
input_stats=$gwas/Table_sumStats_7th_run.csv
Rscript --vanilla $scripts/prepare_summary_stats.R $input_stats $gwas_name

#Read the lead SNPs stats into R and output as various required formats.
input_file=$gwas/Table_signHits_5th_500kb_run.csv
Rscript --vanilla $scripts/prepare_lead_SNPs.R $input_file $gwas_name


