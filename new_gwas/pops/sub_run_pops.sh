#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

pops=/newhome/qh18484/bin/pops

cd $pops

python pops.feature_selection.py \
	--features PoPS.features.txt.gz \
	--gene_results $gwas \
	--out $gwas 

for CHR in {1..22}
do
python pops.predict_scores.py \
	--gene_loc gene_loc.txt \
	--gene_results $gwas \
	--features PoPS.features.txt.gz \
	--selected_features ${gwas}.features\
	--control_features control.features \
	--chromosome $CHR \
	--out $gwas 
done