#!/bin/bash
HOME=/mnt/storage/home/qh18484
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/loci_definition
gwas=$HOME/scratch/new_gwas/gwas/raw
analysis=$HOME/scratch/new_gwas/loci_definition

cd $gwas

gwas_name="eczema21_discovery"

#Read the summary stats into R, rename columns so that consistent with Paternoster2015 files and output as TSV.
input_stats=$gwas/Table_sumStats_7th_run_updated_RSID.csv
Rscript --vanilla $scripts/prepare_summary_stats.R $input_stats $gwas_name

#Filter out variants with no rsid for FUMA
(head -1 results.${gwas_name}_fuma.txt; tail -n +2 results.${gwas_name}_fuma.txt | awk -v OFS="\t" '($7 ~ "rs") {print $0}') > results.${gwas_name}_fuma_ver2.txt 

#Read the lead SNPs stats into R and output as various required formats.
input_file=$gwas/Table_signHits_5th_500kb_run.csv
Rscript --vanilla $scripts/prepare_lead_SNPs.R $input_file $gwas_name

#Basic data check - in each column in summary data, count the number of values and explore the top hits.
Rscript --vanilla $scripts/gwas_data_check.R $gwas_name

cd $analysis

#LiftOVer the GRCH37 intervals from MendelVar to GRCH38, used for colocalization with eQTL Catalogue data.
$HOME/bin/liftover/liftOver ${gwas_name}_interval.bed $HOME/bin/liftover/hg19ToHg38.over.chain \
${gwas_name}_interval_mapped.bed ${gwas_name}_interval_unmapped.bed

#All intervals mapped successfully