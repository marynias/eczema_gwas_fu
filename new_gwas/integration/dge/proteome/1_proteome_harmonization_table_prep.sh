#!/bin/bash
analysis=$HOME/scratch/new_gwas/integration/dge/proteome
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/dge/proteome
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/scratch/new_gwas/genome
master=$HOME/scratch/new_gwas/integration/1_loci_prep
integration=$HOME/scratch/new_gwas/integration/final_integration

gwas_name="eczema21_discovery"

cd $analysis

###########The below code has already been run. It only needs to be run once.
#Make sure that the IDs are converted to HUGO gene symbol.
: '
python $old_scripts/sync_ids.py --tab Elias2017_TableE2.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Elias2017_TableE2.txt

python $old_scripts/sync_ids.py --tab Molin2015_Table3.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Molin2015_Table3.txt

python $old_scripts/sync_ids.py --tab Morelli2021_exd14276-sup-0006-tables1.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Morelli2021_exd14276-sup-0006-tables1.txt

for a in Pavel2019*
do
python $old_scripts/sync_ids.py --tab $a \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp $a
done

Rscript --vanilla $scripts/prepare_proteome_tables.R
'

Rscript --vanilla $scripts/merge_proteome_tables.R  $master/${gwas_name}_master.csv \
${gwas_name}_dge_proteome.csv

#Move final output file to the final integration folder.
cp ${gwas_name}_dge_proteome.csv $integration
