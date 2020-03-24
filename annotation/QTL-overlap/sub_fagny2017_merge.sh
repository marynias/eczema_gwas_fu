#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/annotation/QTL-overlap]
data_manipulation=$HOME/analysis/annotation/data_manipulation

cd $PBS_O_WORKDIR

Rscript --vanilla $scripts/fagny2017_merge.R whole_blood $data_manipulation/interval_r2_0.2_1k.snps
Rscript --vanilla $scripts/fagny2017_merge.R skin $data_manipulation/interval_r2_0.2_1k.snps