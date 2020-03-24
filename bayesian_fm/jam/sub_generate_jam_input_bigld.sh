#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam

cd $PBS_O_WORKDIR

python $scripts/generate_jam_input3.py --tab $locus \
--vcf $my_vcf --out chr${chr}.${snp}.bigld --pos 3 --ref 4 \
--alt 5 --beta 8 --se 9 --rsid 12