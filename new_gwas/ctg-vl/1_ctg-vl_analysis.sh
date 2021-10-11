#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/ctg-vl/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/ctg-vl
gwas=$HOME/data/gwas/paternoster2015

cd $analysis

#Prepare CTG-VL input
#The file should be tab separated and minimum should contain the columns CHR, BP, SNP, A1, A2, FREQ, BETA, SE and P. The column "N" is optional 
cat $gwas/results.euro.pval.all.1k | sed 's/POS/BP/' | sed 's/EFFECT_ALLELE/A1/' | sed 's/NON_EFFECT_ALLELE/A2/' | sed 's/EFFECT_ALLELE_FREQ/FREQ/' | sed 's/PVAL/P/' | awk -v OFS="\t" '{print $2, $3, $1, $4, $5, $7, $8, $9, $11, $6}' >results.euro.pval.all.ctg

#Filter out variants with no rsid and not ATCG.
(head -1 results.euro.pval.all.ctg; tail -n +2 results.euro.pval.all.ctg | awk -v OFS="\t" '($3 ~ "rs") {print $0}') >results.euro.pval.all.ctg_ver2

###For analysis involving genetic correlation, remove the MHC region.
##Filter out the MHC region
##GRCh37
#https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
#chr6:28,477,797-33,448,354

##GRCh38
#https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
#chr6:28,510,120-33,480,577

cat $gwas/results.euro.pval.all.1k | sed 's/POS/BP/' | sed 's/EFFECT_ALLELE/A1/' | sed 's/NON_EFFECT_ALLELE/A2/' | sed 's/EFFECT_ALLELE_FREQ/FREQ/' | sed 's/PVAL/P/' | awk -v OFS="\t" '($2 == 6 && $3 > 3344835 || $3 < 2847779) || ($2 != 6) {print $2, $3, $1, $4, $5, $7, $8, $9, $11, $6}' >results.euro.pval.nomhc.ctg

#Filter out variants with no rsid and not ATCG.
(head -1 results.euro.pval.all.ctg; tail -n +2 results.euro.pval.all.ctg | awk -v OFS="\t" '($3 ~ "rs") {print $0}') >results.euro.pval.nomhc.ctg_ver2

