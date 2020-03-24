#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/annotation/fitCons

module load languages/R-3.5.1-ATLAS-gcc-6.1

cd $PBS_O_WORKDIR

Rscript $scripts/annotate_fitcons.R