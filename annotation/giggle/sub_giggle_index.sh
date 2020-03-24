#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
giggle=$HOME/bin/giggle

cd $PBS_O_WORKDIR

cat $input | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > ${input%.bed}.bed.gz 
$giggle/bin/giggle index -f -i ${input%.bed}.bed.gz -o $output 