#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/depict/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/depict
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015

cd $analysis

#We recommend doing 2 DEPICT runs based on independent genome-wide significant SNPs 
#and another based on all independent SNPs with P value < 1e-5.

#####Top independent loci analysis
#Top p-value 
cat $gwas/paternoster_2015_index_snps_sorted.txt | cut -f1 >paternoster_2015_depict.top

######All independent SNPs with P value < 1e-5
cat $gwas/results.euro.pval.1k.dbsnp |  awk '($11 < 0.00001) {print $12}' >paternoster_2015_depict.1e5

#Before every run change parameters within the depict.py file in /panfs/panasas01/sscm/qh18484/bin/depict
qsub $scripts/sub_run_depict.sh