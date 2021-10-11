#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

cd /panfs/panasas01/sscm/qh18484/bin/depict/

./depict.py