#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
HOME=/panfs/panasas01/sscm/qh18484
cd $HOME/analysis/bayesian_fm/jam
# run your program
Rscript --vanilla $HOME/bin/R2BGLiMS/Examples/JAM_Examples.R
