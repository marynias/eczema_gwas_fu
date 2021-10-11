#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/ldhub/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/ldhub
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015

cd $analysis

#Prepare LDhub input
cat $gwas/results.euro.pval.all.1k | sed 's/EFFECT_ALLELE/A1/' | sed 's/NON_EFFECT_ALLELE/A2/' | sed 's/PVAL/P-value/' | sed 's/SNP/snpid/' | sed 's/Z_SCORE/Zscore/' | awk -v OFS="\t" '{print $1, $4, $5, $10, $6, $11}' >results.euro.pval.all.ldhub.txt

#Filter out variants with no rsid and not ATCG.
(head -1 results.euro.pval.all.ldhub.txt; tail -n +2 results.euro.pval.all.ldhub.txt | awk -v OFS="\t" '($1 ~ "rs") {print $0}') >results.euro.pval.all.ldhub_ver2.txt

zip results.euro.pval.all.ldhub_ver2.txt.zip results.euro.pval.all.ldhub_ver2.txt