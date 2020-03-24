#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
HOME=/panfs/panasas01/sscm/qh18484
FINEMAP_DIR=$HOME/bin/finemap_v1.3_x86_64
cd $HOME/analysis/bayesian_fm/finemap
# run your program

HOME=/panfs/panasas01/sscm/qh18484
FINEMAP_DIR=$HOME/bin/finemap_v1.3_x86_64
cd $FINEMAP_DIR
# run your program

./finemap_v1.3_x86_64 --sss --in-files master --dataset 1


