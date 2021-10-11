#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/garfield/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/garfield
garfield=$HOME/bin/garfield/garfield-v2

#The pipeline requires GChr37 coordinates.

#Preparation of Paternoster 2015 GWAS for input. Full summary stats SNPs.
#Need to use a header-less file with summary stats.

#First modify the input generating file in $HOME/bin/garfield/garfield-v2/garfield-create-input-gwas.sh
$HOME/bin/garfield/garfield-v2/garfield-create-input-gwas.sh

#Now modify the input folder in the garfield run file $HOME/bin/garfield/garfield-v2/garfield
cd $garfield
qsub $scripts/sub_run_garfield.sh
