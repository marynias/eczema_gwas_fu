#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/fitCons
scripts=$HOME/bin/eczema_gwas_fu/annotation/fitCons
prog_dir=$HOME/bin/kggseq10hg19
data_manipulation=$HOME/analysis/annotation/data_manipulation

cd $analysis

#Annotate our SNPs with their fitCons score.
qsub $scripts/sub_fitcons.sh