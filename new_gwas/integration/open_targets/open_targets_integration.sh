#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/open_targets
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/open_targets
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
for my_file in *tsv
do
python $old_scripts/sync_ids.py --tab $my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file%.tsv}.txt
done
#Annotate the master file.
Rscript --vanilla $scripts/open_targets_annotation.R $master/paternoster2015_master.csv \
paternoster2015_open_targets.csv
