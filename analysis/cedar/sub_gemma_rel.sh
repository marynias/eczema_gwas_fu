#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR

gemma -g $my_bimbam -p $my_pheno -gk 1 -o $my_output