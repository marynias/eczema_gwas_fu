#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
sample_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/sample

cd $PBS_O_WORKDIR

plink --gen data.all.gen --sample $sample_data/data.chr01.sample --oxford-single-chr 1 --make-bed --out data.all