#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam

cd $PBS_O_WORKDIR

Rscript $scripts/run_jam2.R chr${chrom}.$snp.$interval.beta chr${chrom}.$snp.$interval.matrix chr${chrom}_${snp}_${interval}_${r2}r2_${maf}maf.results