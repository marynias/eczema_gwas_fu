#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=48:00:00
#PBS -V

cd $PBS_O_WORKDIR

vcftools=/panfs/panasas01/sscm/qh18484/bin/vcftools/bin
$vcftools/vcftools --gzvcf $my_file --geno-r2 --stdout | gzip -c > ${my_file%.vcf.gz}.ld.gz

