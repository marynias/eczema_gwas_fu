#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR

%%	my_file=do_wymiany
%%	table=do_wymiany

to_remove=$(awk '{print "--remove-indv " $1}' $table | tr '\n' ' ')

#Remove all monomorphic SNPs and remove listed individuals
vcftools=/panfs/panasas01/sscm/qh18484/bin/vcftools/bin
$vcftools/vcftools --gzvcf $my_file $to_remove --mac 1 --recode --stdout | gzip -c > ${my_file%.vcf.gz}_no_rels_no_mono.vcf.gz
