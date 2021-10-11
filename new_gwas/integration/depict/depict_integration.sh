#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/depict
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/depict
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
python $old_scripts/sync_ids.py --tab paternoster_2015_depict.top_geneprioritization.txt --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >paternoster_2015_depict.top_geneprioritization_syn.txt
#Annotate the master file.
Rscript --vanilla $scripts/depict_annotation.r $master/paternoster2015_master.csv \
paternoster_2015_depict.top_geneprioritization_syn.txt  \
paternoster2015_depict.csv

