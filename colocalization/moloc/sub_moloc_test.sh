#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00

HOME=/panfs/panasas01/sscm/qh18484
cd $HOME/analysis/colocalization/moloc/data

# run your program
Rscript --vanilla $HOME/bin/eczema_gwas_fu/colocalization/moloc/moloc_test.R