#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/pops
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/pops
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
python $old_scripts/sync_ids.py --tab paternoster2015_pops_annotated.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >paternoster2015_pops_annotated_syn.txt
#Annotate the master file.
Rscript --vanilla $scripts/pops_annotation.R $master/paternoster2015_master.csv \
paternoster2015_pops_annotated_syn.txt  \
paternoster2015_pops.csv
