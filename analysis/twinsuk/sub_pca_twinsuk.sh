#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk

cd $PBS_O_WORKDIR

Rscript --vanilla $scripts/pca_twinsuk.R