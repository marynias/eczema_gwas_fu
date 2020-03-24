#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel

cd $PBS_O_WORKDIR

python $scripts/generate_BigLD_input.py $file ${file%.vcf*}