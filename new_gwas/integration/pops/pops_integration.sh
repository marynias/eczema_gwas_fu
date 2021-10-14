#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/integration/pops
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/pops
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/scratch/new_gwas/genome
master=$HOME/scratch/new_gwas/integration/1_loci_prep
pops=$HOME/scratch/pops/
pops_analysis=$HOME/scratch/new_gwas/pops

gwas_name="eczema21_discovery"

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
python $old_scripts/sync_ids.py --tab $pops_analysis/${gwas_name}_pops_annotated.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${gwas_name}_pops_annotated_syn.txt

#Annotate the master file.
Rscript --vanilla $scripts/pops_annotation.R $master/${gwas_name}_master.csv \
${gwas_name}_pops_annotated_syn.txt  \
${gwas_name}_pops.csv
