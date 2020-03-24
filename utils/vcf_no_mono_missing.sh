#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR

%%	my_file=do_wymiany
%%	maxm=do_wymiany
%%	output=do_wymiany

#Keep only inds listed in the table and remove all monomorphic SNPs
vcftools=/panfs/panasas01/sscm/qh18484/bin/vcftools/bin
$vcftools/vcftools --gzvcf $my_file --max-missing $maxm --mac 1 --recode --stdout | gzip -c > $output