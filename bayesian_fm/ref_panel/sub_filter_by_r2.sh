#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts_ref=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015

cd $PBS_O_WORKDIR

my_position=$(grep -w $my_snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -d$'\t' -f3)
python $scripts_ref/filter_by_r2.py $my_position $my_r2 $my_ld_file >${my_snp}_${my_position}_15000_${my_r2}.${my_dataset}