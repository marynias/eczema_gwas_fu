#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
HOME=/panfs/panasas01/sscm/qh18484
cd $PBS_O_WORKDIR

$HOME/bin/finemap_v1.3_x86_64/finemap_v1.3_x86_64 --sss --in-files master