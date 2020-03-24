#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas
twas=$HOME/bin/fusion_twas-master

PATH=$PATH:$twas
cd $PBS_O_WORKDIR

Rscript $twas/FUSION.compute_weights.R \
--bfile $input \
--tmp $temp \
--out $output \
--models top1,blup,lasso,enet \
--PATH_plink /cm/shared/apps/Plink2/plink \
--PATH_gcta $twas/gcta_nr_robust \
--PATH_gemma /cm/shared/apps/GEMMA/GEMMA-0.96/bin/gemma \
--verbose 2