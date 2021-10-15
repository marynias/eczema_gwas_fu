#!/bin/bash
HOME=/mnt/storage/home/qh18484
snpsea=$HOME/scratch/snpsea
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/snpsea
analysis=$HOME/scratch/new_gwas/snpsea
gwas=$HOME/scratch/new_gwas/gwas/raw
data_manipulation=$HOME/scratch/new_gwas/genome

cd $analysis

gwas_name="eczema21_discovery"

#Prepare SNPSea input
(tail -n +2 $gwas/leadSNPs.${gwas_name}.txt | cut -f2,3,12,11 | awk -v OFS="\t" '{print $1, $2, $4, $3}' | 
sed 's/^/chr/' | sort -k1,1V -k2,2n ) | sed '1 i\CHR\tPOS\tSNP\tP' >${gwas_name}.snpsea

sbatch --export=ALL,my_gwas=${gwas_name}.snpsea,gene_matrix=$snpsea/GeneAtlas2004.gct.gz,output=${gwas_name}_GeneAtlas2004 $scripts/sub_snpsea.sh
sbatch --export=ALL,my_gwas=${gwas_name}.snpsea,gene_matrix=$snpsea/ImmGen2012.gct.gz,output=${gwas_name}_ImmGen2012 $scripts/sub_snpsea.sh
sbatch --export=ALL,my_gwas=${gwas_name}.snpsea,gene_matrix=$snpsea/FANTOM2014.gct.gz,output=${gwas_name}_FANTOM2014 $scripts/sub_snpsea.sh
sbatch --export=ALL,my_gwas=${gwas_name}.snpsea,gene_matrix=$snpsea/GO2013.gct.gz,output=${gwas_name}_GO2013 $scripts/sub_snpsea.sh

source activate myenv

#Visualize the results
for a in ${gwas_name}_GeneAtlas2004 ${gwas_name}_ImmGen2012 ${gwas_name}_FANTOM2014 ${gwas_name}_GO2013
do
sbatch --export=ALL,output=$a $scripts/sub_snpsea_viz.sh
done

#Sort the results
for a in ${gwas_name}_GeneAtlas2004 ${gwas_name}_ImmGen2012 ${gwas_name}_FANTOM2014 ${gwas_name}_GO2013
do
(head -1 $a/condition_pvalues.txt; tail -n +2 $a/condition_pvalues.txt | sort -t$'\t' -k2n) > $a/condition_pvalues_sorted.txt 
done