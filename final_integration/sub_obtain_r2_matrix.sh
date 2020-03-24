#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=04:00:00
#PBS -V

cd $PBS_O_WORKDIR


HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/final_integration

python $scripts/obtain_r2_matrix.py ${my_snp}_snps.unadjusted_Model09.interval.gassocplot \
$HOME/working/data/results/RefPanel/1kGenomes/${my_snp}_3Mbp_1kEUR.ld.gz ${my_snp}_gassocplot.ld