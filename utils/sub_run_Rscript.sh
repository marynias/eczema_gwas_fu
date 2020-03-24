#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=120:00:00
#PBS -V

cd $PBS_O_WORKDIR

Rscript --vanilla $my_script
