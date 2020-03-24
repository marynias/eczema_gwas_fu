#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
snpsea=$HOME/bin/snpsea

cd $PBS_O_WORKDIR

python $snpsea/bin/snpsea-barplot $output
python $snpsea/bin/snpsea-heatmap $output
Rscript $snpsea/bin/snpsea-type1error $output