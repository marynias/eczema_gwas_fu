#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/smultixcan
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/smultixcan
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

for my_file in *csv
do
python $old_scripts/sync_ids.py --tab $my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 2 --delim $',' --delim_c $',' --convert upper >${my_file%.txt}.txt
done

#Annotate the master file.
Rscript --vanilla $scripts/smultixcan_annotation.R $master/paternoster2015_master.csv \
paternoster2015_smultixcan.csv
