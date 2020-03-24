#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/annotation/QTL-overlap

cd $PBS_O_WORKDIR

Rscript $scripts/fagny2017_extract.R whole_blood
Rscript --vanilla $scripts/fagny2017_extract.R skin