#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scrach/new_gwas/ldhub/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/ldhub
gwas=$HOME/scratch/new_gwas/gwas/raw

cd $analysis

gwas_name="eczema21_discovery"

#Prepare LDhub input
cat $gwas/results.${gwas_name}.txt| sed 's/EFFECT_ALLELE/A1/' | sed 's/NON_EFFECT_ALLELE/A2/' | sed 's/PVAL/P-value/' | sed 's/RSID/snpid/' | sed 's/Z_SCORE/Zscore/' | awk -v OFS="\t" '{print $12, $4, $5, $10, $6, $11}' >results.${gwas_name}.ldhub.txt

#Filter out variants with no rsid and not ATCG.
(head -1 results.${gwas_name}.ldhub.txt; tail -n +2 results.${gwas_name}.ldhub.txt | awk -v OFS="\t" '($1 ~ "rs") {print $0}') >results.${gwas_name}.ldhub_ver2.txt

zip results.${gwas_name}.ldhub_ver2.txt.zip results.${gwas_name}.ldhub_ver2.txt