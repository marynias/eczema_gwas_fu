#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc

cd $PBS_O_WORKDIR

python $scripts/extract_gtex.py $gtex_data $ensembl_selection $gene_names $gtex_ref $tissue