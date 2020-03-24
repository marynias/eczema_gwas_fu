#!/bin/bash
#set -eu
#set -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/final_integration/gene
scripts=$HOME/bin/eczema_gwas_fu/final_integration
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
datasets=$HOME/working/data/Datasets
results=$HOME/working/data/results/euro

cd $analysis

#Collect results from coloc and TWAS.
#Coloc genes and SNPs.
for a in $results/coloc/raw/blueprint
do
for var in mono tcel neut
do
my_file=${var}.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 6 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Blueprint_${var}_coloc" --study_id_column "Blueprint" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 5 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 197 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Blueprint_${var}_coloc" --study_id_column "Blueprint" --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 6 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 7 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --sample_size 197 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done

for a in $results/coloc/raw/cedar
do
for var in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
my_file=cedar_${var}.coloc_results.b
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 6 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "CEDAR_${var}_coloc" --study_id_column CEDAR --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 5 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 323 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "CEDAR_${var}_coloc" --study_id_column CEDAR --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 6 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 7 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --sample_size 323 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done

for a in $results/coloc/raw/eqtlgen
do
my_file=eqtlgen.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "eQTLgen_eQTL_coloc" --study_id_column "eQTLgen" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 4 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue "whole blood" --sample_size 14115 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "eQTLgen_eQTL_coloc" --study_id_column "eQTLgen" --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 5 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 6 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "whole blood" --sample_size 14115 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done


for a in $results/coloc/raw/gtex
do
for var in Skin_Not_Sun_Exposed_Suprapubic \
Spleen \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes
do
my_file=${var}.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 6 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "GTEx_${var}_coloc" --study_id_column "GTEx" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 4 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 300 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "GTEx_${var}_coloc" --study_id_column "GTEx" --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 6 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 7 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --sample_size 300 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done

for a in $results/coloc/raw/kim-hellmuth
do
for var in lps6h lps90 mdp6h mdp90 ctrl
do
my_file=kh_${var}.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Kim-Hellmuth_${var}_coloc" --study_id_column "Kim-Hellmuth" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 4 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 185 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Kim-Hellmuth_${var}_coloc" --study_id_column "Kim-Hellmuth"  --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 5 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 6 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --sample_size 185 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done


for a in $results/coloc/raw/sun
do
my_file=sun.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Sun_pQTL_coloc" --study_id_column "Sun" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 4 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue "whole blood" --sample_size 3301 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "Sun_pQTL_coloc" --study_id_column "Sun" --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 5 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 6 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "whole blood" --sample_size 3301 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done


for a in $results/coloc/raw/twinsuk
do
for var in lcl skin 
do
for sample in eczema noneczema all 
do
my_file=${sample}_${var}.coloc_results
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 6 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "TwinsUK_${sample}_${var}_coloc" --study_id_column "TwinsUK" --output ".gene" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 5 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 700 --type "coloc - gene" --cis_trans cis --evidence_weight 1  --enable_counts 1

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "TwinsUK_${sample}_${var}_coloc" --study_id_column "TwinsUK" --output ".snp" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 6 --gene 0 --index_locus 2 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 7 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --sample_size 700 --type "coloc - SNP" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done 
done

for a in $results/twas/cedar
do
for var in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
my_file=cedar_${var}.all

#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string "CEDAR_${var}_TWAS" --study_id_column "CEDAR" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 24 --index_locus_string "" --beta 0 --beta_se 20 --posterior_prob 5 --bayes_factor 0 --score 0 \
--pvalue 19 --fdr 23 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 323 --type "TWAS" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done



for a in $results/twas/cedar
do
for var in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
my_file=cedar_${var}.coloc.sig0.45

#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string "CEDAR_${var}_TWAS_coloc"  --study_id_column "CEDAR" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 26 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 25 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 323 --type "TWAS coloc" --cis_trans cis --evidence_weight 1  --enable_counts 1
done
done

my_file=gtex_twas.all
python $scripts/sync_ids.py --tab $results/twas/gtex/${my_file} --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string 1 --study_id_column "GTEx" --output "GTEx_TWAS" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 24 --index_locus_string "" --beta 19 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 23 --fdr 20 --sig_threshold 0.05 --effect_allele 0 --tissue "PANEL" --sample_size 300 --type "TWAS" --cis_trans cis --evidence_weight 1  --enable_counts 1

sed -i 's/^/GTEx_TWAS_/' GTEx_TWAS.processed

my_file=coloc_twas.sig0.45
python $scripts/sync_ids.py --tab $results/twas/gtex/${my_file} --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string 1 --study_id_column "GTEx" --output "GTEx_TWAS_coloc"  --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 26 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 25 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue "PANEL" --sample_size 300 --type "TWAS coloc" --cis_trans cis --evidence_weight 1  --enable_counts 1

sed -i 's/^/GTEx_TWAS_coloc_/' GTEx_TWAS_coloc.processed

for a in $results/twas/twinsuk
do
for var in skin lcl
do
my_file=twinsuk_${var}.all

#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string 1 --study_id_column "TwinsUK" --output "TwinsUK_${var}_TWAS" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 24 --index_locus_string "" --beta 20 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 19 --fdr 23 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 700 --type "TWAS" --cis_trans cis --evidence_weight 1  --enable_counts 1
sed -i "s/^/TwinsUK_${var}_TWAS_/" TwinsUK_${var}_TWAS.processed
done
done

for a in $results/twas/twinsuk
do
for var in skin lcl
do
my_file=twinsuk_${var}.coloc.sig0.45

#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string 1 --study_id_column "TwinsUK" --output "TwinsUK_${var}_TWAS_coloc" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 26 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 25 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.45 --effect_allele 0 --tissue $var --sample_size 700 --type "TWAS coloc" --cis_trans cis --evidence_weight 1  --enable_counts 1
sed -i "s/^/TwinsUK_${var}_TWAS_coloc_/" TwinsUK_${var}_TWAS_coloc.processed
done
done

a=$HOME/analysis/annotation/networks/prixfixe
my_file=non_redundant_loci.prixfixe.out

python $scripts/sync_ids.py  --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "PrixFixe" --study_id_column "PrixFixe" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 9 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --sample_size "" --type "network analysis" --cis_trans "" --evidence_weight 3  --enable_counts 0


a=$datasets/eczema/DGE/Cole2014/Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp1.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 8 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 36 --type "DGE" --cis_trans "" --evidence_weight 2 --enable_counts 0

a=$datasets/eczema/DGE/Cole2014
for my_file in Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp3.txt \
Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp4.txt \
Cole_\(2014\)_Filaggrin-stratified_transcriptomic_analysis_supp5.txt 
do
python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a/$my_file \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 9 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 36 --type "DGE" --cis_trans "" --evidence_weight 2 --enable_counts 0
done

a=$datasets/eczema/DGE/Ewald2015/12920_2015_133_MOESM7_ESM.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 0 --beta 5 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 7 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 97 --type "DGE meta" --cis_trans "" --evidence_weight 2 --enable_counts 0

a=$datasets/eczema/DGE/Ewald2015/ewald_2015_tab2_MTGDR_ad_classifier.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 2 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 97 --type "MTGDR classifier" --cis_trans "" --evidence_weight 2 --enable_counts 0


a=$datasets/eczema/DGE/Saaf2008/Saaf2008_Table_S3.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 4 --index_locus 0 --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 11 --type "DGE" --cis_trans "" --evidence_weight 2 --enable_counts 0


a=$datasets/eczema/DGE/Ghosh2015
for my_file in Ghosh2015_S1_upregulated_genes.txt \
Ghosh2015_S2_downregulated_genes.txt
do
python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a/$my_file \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 2 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 3 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size "" --type "DGE meta" --cis_trans "" --evidence_weight 2 --enable_counts 0
done


a=$datasets/eczema/protein/Elias2017/Elias2017_TableE2.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 0 --beta 4 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 5 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 10 --type "proteome" --cis_trans "" --evidence_weight 2 --enable_counts 0

a=$datasets/eczema/protein/Molin2015/Molin2015_Table3.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 4 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 5 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "hand skin" --sample_size 12 --type "proteome" --cis_trans "" --evidence_weight 2 --enable_counts 0

a=$datasets/eczema/methylation/Rodriguez2014
for my_file in Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-Supp_Tabl1_AL_vs_AN.txt \
Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-Supp_Tabl1_AL_vs_NN.txt \
Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-Supp_Tabl1_AN_vs_NN.txt 
do
python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a/$my_file \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 4 --index_locus 0 --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 57 --type "methylation" --cis_trans "" --evidence_weight 2 --enable_counts 0
done

for my_file in Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-pages-31-55_Supp_Tab4_AL_vs_AN.txt \
Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-pages-31-55_Supp_Tab4_AL_vs_NN.txt \
Rodriguez_\(2014\)_An_Integrated_Epigenetic_and_Transcriptomic_Analysis_supp-pages-31-55_Supp_Tab4_AN_vs_NN.txt
do
python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a/$my_file \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 5 --index_locus 0 --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 57 --type "DGE" --cis_trans "" --evidence_weight 2 --enable_counts 0
done

a=$datasets/eczema/methylation/Quraishi2015/Quraishi_\(2015\)_Identifying_CpG_sites_supp_sig.txt
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 4 --index_locus 0 --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 366 --type "methylation" --cis_trans "" --evidence_weight 2 --enable_counts 0

a=$datasets/eczema/DGE/Winge2011
for my_file in Supp_Tab_S1_flg.txt \
Supp_Tab_S1_flg_het.txt \
Supp_Tab_S1_flg_no.txt 
do
python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 1 --index_locus 0 --beta 4 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 3 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "skin" --sample_size 58 --type "DGE" --cis_trans "" --evidence_weight 2 --enable_counts 0
done

#regfm results
for a in $results/annotation/regfm/paternoster2015_1k/r2_min0.1/Final-Results.txt
do
my_file=$(basename $a)
#Filter out only significant results and tissues of no interest.
#(head -1 $a; awk -F "\t" -v OFS="\t" '$18 < 0.05' $a) > ${my_file%.txt}.filtered
python $scripts/make_input_table_ver2.py --tab ${my_file%.txt}.filtered  \
--tab_id_string "regfm_r2_min0.1" --study_id_column regfm  --output "regfm_r2_min0.1" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 17 --index_locus 5 --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 18 --sig_threshold 0.05 --effect_allele 0 --tissue Cell  --type "regfm" --cis_trans "" --evidence_weight 1 --enable_counts 1
done

for a in $datasets/hQTL/Pelikan2018/41467_2018_5328_MOESM6_ESM.1.tsv
do
my_file=$(basename $a)
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 4 --index_locus 0  --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 5 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "LCL" --sample_size 25 --type "hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 0
done

for a in $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM3_ESM
do
my_file=$(basename $a)
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 6 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 6 --index_locus 0  --index_locus_string "" --beta 14 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 15 --sig_threshold 0 --effect_allele 0 --tissue "activated T cells" --sample_size 20 --type "DGE" --cis_trans "" --evidence_weight 2  --enable_counts 0
done

for a in $datasets/promoter-enhancer/Rubin2017/ng.3935-S6.1_gained.tsv \
$datasets/promoter-enhancer/Rubin2017/ng.3935-S6.1_stable.tsv \
$datasets/promoter-enhancer/Rubin2017/ng.3935-S6.1_stable_gained.tsv
do
my_file=$(basename $a)
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 0 --my_id 1 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 0  --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "keratinocytes" --sample_size "" --type "promoter enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 0
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc2.1.tsv
do
my_file=$(basename $a)
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 3 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 0 --gene 3 --index_locus 0  --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 12 --fdr 13 --sig_threshold 0 --effect_allele 0 --tissue "keratinocytes" --sample_size "" --type "DEG" --cis_trans "" --evidence_weight 2  --enable_counts 0
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.3_active_epsc.txt
do
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 0 --my_id 1 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".EpiSC_active" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 0  --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC_active" --sample_size "" --type "enhancers" --cis_trans "" --evidence_weight 2  --enable_counts 0
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.3_denovo_keratinocytes.txt
do
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 0 --my_id 1 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".Keratinocytes_denovo" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 0  --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_denovo" --sample_size "" --type "enhancers" --cis_trans "" --evidence_weight 2  --enable_counts 0
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.3_maintained_keratinocytes.txt
do
python $scripts/sync_ids.py --tab $a --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 0 --my_id 1 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".Keratinocytes_maintained" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 0 --gene 1 --index_locus 0  --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_maintained" --sample_size "" --type "enhancers" --cis_trans "" --evidence_weight 2  --enable_counts 0
done
done





