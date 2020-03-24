#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00

HOME=/panfs/panasas01/sscm/qh18484
PAINTOR_DIR=$HOME/bin/PAINTOR_V3.0
cd $HOME/analysis/bayesian_fm/paintor
mkdir SampleData

# run your program
$PAINTOR_DIR/PAINTOR -input $PAINTOR_DIR/SampleData/input.files -in $PAINTOR_DIR/SampleData/ \
-out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS
