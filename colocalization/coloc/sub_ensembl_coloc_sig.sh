#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc

#Annotate results of colocalisation analysis with Ensembl info.
cd $PBS_O_WORKDIR

for a in *.sig*.annotated
do
Rscript $scripts/annotate_ensembl.R $a
done
