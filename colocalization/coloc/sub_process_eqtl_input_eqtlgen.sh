#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015

cd $PBS_O_WORKDIR

while read line
do
set -- $line
gene=$6
echo $line >>${a%.eqtl}_${gene}.eqtl
done <$a


