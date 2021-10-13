#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/depict
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/depict
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/scratch/new_gwas/genome
master=$HOME/scratch/new_gwas/integration/1_loci_prep
depict=$HOME/scratch/depict/results

gwas_name="eczema21_discovery"

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
python $old_scripts/sync_ids.py --tab $depict/${gwas_name}_geneprioritization.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${gwas_name}_geneprioritization_syn.txt
#Annotate the master file.
Rscript --vanilla $scripts/depict_annotation.r $master/${gwas_name}_master.csv \
${gwas_name}_geneprioritization_syn.txt \
${gwas_name}_depict.csv

