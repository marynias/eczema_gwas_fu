#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/moloc
analysis=$HOME/analysis/colocalization/moloc/blood
eqtlgen=$HOME/analysis/colocalization/coloc/eqtlgen/euro/gene
sun=$HOME/analysis/colocalization/coloc/sun_pqtl

#Moloc analysis on cis-eQTLs from EQTLgen and pQTLs from Sun 2018.
cd $analysis

for a in $sun/*.sun
do
my_filename=$(basename $a)
my_gwas=$(echo $my_filename | cut -d"_" -f1,2)
Rscript --vanilla $scripts/moloc_eqtlgen_sun2018.R $sun/${my_gwas}.euro.pval $eqtlgen/${my_filename%.sun}.eqtl $a 
done

#Find out if any high values of colocalisation found. 
for a in *.pp.moloc
do
awk 'NR == 15 && $5 > 0.05 {print FILENAME}' $a
done
