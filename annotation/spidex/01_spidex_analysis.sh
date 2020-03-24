#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/spidex
scripts=$HOME/bin/eczema_gwas_fu/annotation/spidex
utils=$HOME/bin/eczema_gwas_fu/utils
prog_dir=$HOME/bin/kggseq10hg19
data_manipulation=$HOME/analysis/annotation/data_manipulation
spidex=$HOME/working/data/Datasets/splicing/Xiong2014
gwas=$HOME/data/gwas/paternoster2015

cd $analysis

#Pull out significant splicing disruption score (>5) and abs(z-score) > 2. As in http://tools.genes.toronto.edu/
cat $spidex/hg19_spidex.txt | awk -v OFS="\t" '$6 > 5 && ($7 > 2 || $7 < -2) {print}' >$spidex/hg19_spidex_sig.txt
#Annotate variants with 1k EUR-based rsids.

#Look for overlaps between our interval SNPs.
python $utils/update_rsid.py --bim $gwas/results.euro.1k.bim --tab $spidex/hg19_spidex_sig.txt \
--head Y --chrom 1 --pos 2 --ref 4 --alt 5 >$spidex/hg19_spidex_sig_rsid.txt