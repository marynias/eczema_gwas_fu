#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -V

cd /panfs/panasas01/sscm/qh18484/bin/regfm
./regfm.sh $input