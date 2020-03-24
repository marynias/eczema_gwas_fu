#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -V

cd $PBS_O_WORKDIR

gene=$(echo $file | cut -d"_" -f1)
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
group=$(echo $file | cut -d"_" -f6 | cut -d"." -f1)
gemma -g ${snp}_${interval}_${group}.bimbam -p $file -n $factors -a ${snp}_${interval}_${group}.snps -c gemma_${group}_covariate.txt -lmm 4 -o ${gene}_${snp}_${interval}_${group}.gemma -k output/${group}_relatedness.cXX.txt