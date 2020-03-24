#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/twinsuk
scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk

cd $PBS_O_WORKDIR

Rscript --vanilla $scripts/peer.R $expression $covariate