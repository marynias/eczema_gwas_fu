#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

cd $PBS_O_WORKDIR

while read line
do
set -- $line
gene=$7
echo $line >>${my_file%.eqtl}_${gene}.eqtl
done <$my_file