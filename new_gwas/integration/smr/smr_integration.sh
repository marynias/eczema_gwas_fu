#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/smr
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/smr
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

for my_file in *txt
do
python $old_scripts/sync_ids.py --tab $my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file%.txt}.tsv
done

#Annotate the master file.
Rscript --vanilla $scripts/smr_annotation.R $master/paternoster2015_master.csv \
paternoster2015_smr.csv
