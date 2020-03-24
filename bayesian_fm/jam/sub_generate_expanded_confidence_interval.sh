#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

cd $PBS_O_WORKDIR

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam

python $scripts/generate_expanded_confidence_interval.py \
$my_ld \
$my_results \
$my_set \
$my_r2 \
$my_output