#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/integration/dge/gxp
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/integration/dge/gxp
old_scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
master=$HOME/new_gwas/integration/1_loci_prep

cd $analysis

#Make sure that the IDs are converted to HUGO gene symbol.
python $old_scripts/sync_ids.py --tab Table_E1_Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp1.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Table_E1_Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp1.txt

python $old_scripts/sync_ids.py --tab Ewald2015_12920_2015_133_MOESM7_ESM.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Ewald2015_12920_2015_133_MOESM7_ESM.txt

for my_file in Winge2011_Supp_Tab_S1_flg.txt \
Winge2011_Supp_Tab_S1_flg_het.txt \
Winge2011_Supp_Tab_S1_flg_no.txt 
do
python $old_scripts/sync_ids.py --tab $my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp
mv ${my_file}.temp $my_file
done

python $old_scripts/sync_ids.py --tab Tsoi2020_acute_chronic_AD_lesions.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Tsoi2020_acute_chronic_AD_lesions.txt

python $old_scripts/sync_ids.py --tab Rojahn2020_TableE9_control_blisters_downregulated.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Rojahn2020_TableE9_control_blisters_downregulated.txt

python $old_scripts/sync_ids.py --tab Rojahn2020_TableE8_AD_control_blisters_upregulated.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Rojahn2020_TableE8_AD_control_blisters_upregulated.txt

python $old_scripts/sync_ids.py --tab Pavel2020_Table_S3_nonlesional_vs_normal.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Pavel2020_Table_S3_nonlesional_vs_normal.txt

python $old_scripts/sync_ids.py --tab Pavel2020_Table_S3_lesional_vs_normal.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Pavel2020_Table_S3_lesional_vs_normal.txt

python $old_scripts/sync_ids.py --tab He2021_lesional_AD_ver_normal.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp He2021_lesional_AD_ver_normal.txt

python $old_scripts/sync_ids.py --tab He2021_nonlesional_AD_ver_normal.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp He2021_nonlesional_AD_ver_normal.txt

python $old_scripts/sync_ids.py --tab Dyjack2018_TableE4_ADlesional_ver_healthy_control.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Dyjack2018_TableE4_ADlesional_ver_healthy_control.txt

python $old_scripts/sync_ids.py --tab Pavel2019_AD_nonlesional_normal_rnaseq.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Pavel2019_AD_nonlesional_normal_rnaseq.txt

python $old_scripts/sync_ids.py --tab Pavel2019_AD_lesional_normal_rnaseq.txt \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >temp
mv temp Pavel2019_AD_lesional_normal_rnaseq.txt

Rscript --vanilla $scripts/prepare_gxp_tables.R
Rscript --vanilla $scripts/merge_gxp_tables.R $master/paternoster2015_master.csv \
paternoster2015_dge_gxp.csv