#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue

cd $PBS_O_WORKDIR


Rscript --vanilla $scripts/eQTL_Catalogue_coloc_susie_local.R $my_rsid
