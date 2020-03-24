#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts_ref=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel

cd $PBS_O_WORKDIR

Rscript --vanilla $scripts_ref/extract_snp_interval_euro.R $reference $snp $interval