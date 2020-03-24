#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00


# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR

gemma -g all.bimbam -p all.pheno -gk 1 -o all_relatedness