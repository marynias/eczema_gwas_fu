#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/open_targets/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/open_targets_genetics
gwas=$HOME/scratch/new_gwas/gwas/raw

gwas_name="eczema21_discovery"

cd $analysis

############ STOPGAP
############
#Preparation of Paternoster 2015 GWAS for input for STOPGAP. Only lead SNPs.
tail -n +2 $gwas/results.${gwas_name}.txt | (echo -e "chromosome\tbase_pair_location\tvariant_id\teffect_allele\tother_allele\tbeta\tstandard_error\tp-value" && awk -v OFS="\t" '{print $2, $3, $12, $4, $5, $8, $9, $11}') >${gwas_name}_input.tsv

echo -e "chromosome\tbase_pair_location\tvariant_id\teffect_allele\tother_allele\tbeta\tstandard_error\tp-value" >${gwas_name}_open_targets_select.tsv 
while read line
do
set -- $line
my_rsid=$3
grep -w $my_rsid ${gwas_name}_input.tsv >>${gwas_name}_open_targets_select.tsv 
done < $gwas/leadSNPs.${gwas_name}.index

#POSTGAP pipeline has some issues with ENSEMBL calls for certain variants.
#Divide the file into parts to find out which variant(s) is at fault.
tail -n +2 ${gwas_name}_open_targets_select.tsv >${gwas_name}_open_targets_headerless.tsv
split -l 10 ${gwas_name}_open_targets_headerless.tsv ${gwas_name}_open_targets_part

#Affix the header back
HEADER=$(head -1 ${gwas_name}_open_targets_select.tsv)
for i in ${gwas_name}_open_targets_part*
do
    sed -i -e "1i$HEADER" "$i"
done

#Run POST GAP using Ubuntu Virtual Machine on Windows. 
for i in ${gwas_name}_open_targets_part*
do
    echo $i
    python POSTGAP.py --database_dir databases_dir --summary_stats $i --disease eczema --output ${i}_postgap.txt --bayesian
done

#The POSTGAP pipeline works fine when input files divided into small chunks.
#However, check if 144 columns in each file as there have been cases of incorrect output.
#Merge the output files.
Rscript --vanilla $scripts/merge_postgap.R $gwas_name



############ OPEN TARGETS V2G PIPELINE
############
###This V2G pipeline only requires a list of lead RSIDs.
Rscript --vanilla $scripts/run_v2g.R $gwas/leadSNPs.${gwas_name}.index $gwas_name

