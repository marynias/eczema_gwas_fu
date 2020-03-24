#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -V

utils=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/utils

cd $PBS_O_WORKDIR


python $utils/update_rsid.py --bim $my_bim --tab $my_tab \
--head $header --chrom $chrom --pos $pos --ref $ref --alt $alt >$out
