#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V


cd $PBS_O_WORKDIR

$HOME/bin/garfield/garfield-v2/garfield