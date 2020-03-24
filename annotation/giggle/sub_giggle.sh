#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
giggle=$HOME/bin/giggle

cd $PBS_O_WORKDIR

$giggle/bin/giggle search -i $input -q $query -v -o | gzip > ${query}_vs_${input}.gz