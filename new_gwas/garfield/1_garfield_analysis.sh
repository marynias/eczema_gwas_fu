#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/garfield/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/garfield
garfield=$HOME/scratch/garfield/garfield-v2
gwas=$HOME/scratch/new_gwas/gwas/raw

gwas_name="eczema21_discovery"

cd $analysis

#The pipeline requires GChr37 coordinates.

#Preparation of GWAS for input. Full summary stats SNPs.
#Need to use a header-less file with summary stats.
tail -n +2 $gwas/results.${gwas_name}.txt >$gwas/results.${gwas_name}_headerless.txt

#First modify the input generating file and run it.
$garfield/garfield-create-input-gwas.sh

#Now modify the input folder in the garfield run file (garfield-v2/garfield)
cd $garfield
sbatch $scripts/sub_run_garfield.sh
