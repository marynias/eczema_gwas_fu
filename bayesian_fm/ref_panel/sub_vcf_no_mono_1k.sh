#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR

table=integrated_call_samples_v3.20130502.ALL.panel
to_remove=$(awk '($3 != "EUR") {print "--remove-indv " $1}' $table | tr '\n' ' ')

#Remove all monomorphic SNPs and keep inds only in EUR population
vcftools=/panfs/panasas01/sscm/qh18484/bin/vcftools/bin
my_file=ALL.chr${PBS_ARRAYID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
$vcftools/vcftools --gzvcf $my_file $to_remove --mac 1 --recode  --stdout | gzip -c > ${my_file%.vcf.gz}_EUR_no_mono.vcf.gz