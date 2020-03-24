#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/final_integration/lookups
scripts=$HOME/bin/eczema_gwas_fu/final_integration
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
datasets=$HOME/working/data/Datasets
giggle_analysis=$HOME/analysis/annotation/giggle
results=$HOME/working/data/results/euro
utils=$HOME/bin/eczema_gwas_fu/utils

cd $analysis
#Prepare the dbSNP alias table-compliant list of SNPs in the r2 interval of the eczema gwas. Thankfully, no index SNP rsid changed.
#Bed like-file with the last column being the index SNP.
python $data_manipulation_scripts/snp_json_substitute.py --tab $data_manipulation/interval_r2_0.2_1k_nodups.map \
--db $data_manipulation/rsid_synonyms.txt --head N  --rsid 4 --delim A >interval_r2_0.2_1k_nodups_sync.map

#Create an annotation of all 3-Mbp intervals around the index SNPs with all the genes within that interval.
Rscript --vanilla $scripts/annotation_rtracklayer.R

#Substitute the ENSG IDs to HUGO gene names, when possible.
python $data_manipulation_scripts/snp_json_substitute.py --tab \
paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper --head Y  --rsid 5 --delim A >paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.hugo

python $data_manipulation_scripts/snp_json_substitute.py --tab \
paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches_ref \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper --head Y  --rsid 5 --delim A >paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches_ref.hugo


cat $datasets/Resource/Resource_index.tsv | sed 's/^\./\/panfs\/panasas01\/sscm\/qh18484\/working\/data\/Datasets/' >$datasets/Resource/Resource_index.txt 

#Script to substitute rsid/ENSG gene names to those in the alias table/Hugo db. Gene names = uppercase, rsid = lowercase
for a in $datasets/eQTL/Wijst2018/41588_2018_89_MOESM6_ESM.2.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 2 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp


#Script to carry out lookups. Attaches one final column - index SNP to whose interval the given SNP belongs.
#rior to this stage, may want to further extract results significant only for a given tissue etc., if table 
#combines results from many experiments.
#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup


#Script to prepare standardised input table.
python $scripts/make_input_table.py --tab ${my_file}.lookup  --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 13 --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 6 --tissue PBMC --sample_size 45 --type eQTL --cis_trans coexpression --evidence_weight 2
done


for a in $datasets/eQTL/Zhernakova2017/exon_level_eQTLs_independent_effects.txt 
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 11 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 11 --index_locus 14 --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue whole_blood --sample_size 2116 --type eQTL --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Zhernakova2017/polyA-ratio_level_eQTLs_independent_effects.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 13 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 13 --index_locus 14 --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue whole_blood --sample_size 2116 --type "polyA ratio level eQTL" --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Zhernakova2017/gene_level_eQTLs_independent_effects_interactions.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 5 --index_locus 14 --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue whole_blood --sample_size 2116 --type "context specific eQTL" --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Zhernakova2017/gene_level_eQTLs.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 5 --index_locus 23 --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 22 --sig_threshold 0.05 --effect_allele 10 --tissue whole_blood --sample_size 2116 --type "eQTL" --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Zhernakova2017/exon-ratio_level_eQTLs_independent_effects.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 13 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 13 --index_locus 14 --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue whole_blood --sample_size 2116 --type "exon-ratio eQTL" --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Fairfax2014/1246949stableS2.2.tsv
do
for f in LPS2 LPS24 Naive IFN
do
my_file=$(basename $a)
python $scripts/sync_ids.py --tab 1246949stableS2.2_${f}.txt --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper > "1246949stableS2.2_${f}.temp"

python $scripts/do_lookups.py --tab "1246949stableS2.2_${f}.temp" --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 1 --delim $'\t' --delim_c $'' >"1246949stableS2.2_${f}.lookup"

python $scripts/make_input_table.py --tab "1246949stableS2.2_${f}.lookup" --input_path $a \
--tab_id Resource_index.txt --output "_${f}" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 11 --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 9 --sig_threshold 0.05 --effect_allele 0 --tissue monocytes+${f} --sample_size 432 --type eQTL --cis_trans cis --evidence_weight 2
done
done

#Prepare input from Schmiedel 2018
for f in B_CELL_NAIVE CD4_NAIVE CD4_STIM CD8_NAIVE CD8_STIM M2 MONOCYTES NK TFH TH17 TH1 TH2 THSTAR TREG_MEM TREG_NAIVE
do
#cat $datasets/eQTL/Schmiedel2018/${f}.vcf | grep -v "^#" | sed 's/;/\t/g' | sed 's/Gene=//' | sed 's/GeneSymbol=//' | sed 's/Pvalue=//' | sed 's/Beta=//'>$datasets/eQTL/Schmiedel2018/${f}.tsv
#python $scripts/sync_ids.py --tab $datasets/eQTL/Schmiedel2018/${f}.tsv --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $',' --convert upper > ${f}.temp

#python $scripts/do_lookups.py --tab ${f}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 3 --delim $'\t' --delim_c $'' >${f}.lookup

python $scripts/make_input_table.py --tab ${f}.lookup --input_path $datasets/eQTL/Schmiedel2018/${f}.vcf \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 3 --gene 8 --index_locus 12 --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 10 --fdr 0 --sig_threshold 0.05 --effect_allele 5 --tissue ${f} --sample_size 92 --type eQTL --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Pala2017/isoQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 14 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 13 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 14 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 14 --gene 13 --index_locus 15 --beta 9 --beta_se 10 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue whole_blood --sample_size 606 --type eQTL --cis_trans cis --evidence_weight 2
done

for a in $datasets/eQTL/Raj2014/tableS4_eu_cd4T_cis_fdr05.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 11 --beta 9 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 10 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "CD4+ T-cells" --sample_size 401 --type eQTL --cis_trans cis --evidence_weight 2
done

cp *.processed ./to_analyze
cd to_analyze 
Rscript --vanilla $scripts/lookups_integration.R
python $data_manipulation_scripts/snp_json_substitute.py --tab paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.secondary \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper --head Y  --rsid 5 --delim A >paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.secondary.hugo

#Giggle range look-ups.
for a in $datasets/promoter-enhancer/Gasperini2019/mmc2.2-3_combined.txt
do
my_file=$(basename $a)
input_file=${my_file%.txt}.bed
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${input_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${input_file}_index_processed.bed
#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 9 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 9 --index_locus 17 --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "K562 erythroleukemia cell line" --type enhancer --cis_trans "" --evidence_weight 2
done

for a in $datasets/promoter-enhancer/Mumbach2017/ng.3963-S5.1.tsv 
do
my_file=$(basename $a)
for value in x y
do
proc=${my_file%.tsv}_${value}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#sed -i -r 's/ \(\+[0-9]+\)//g' $processed_giggle_output 
#sed -i -r 's/ \(\-[0-9]+\)//g' $processed_giggle_output 

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table.py --tab ${proc}.lookup --input_path $a \
--tab_id Resource_index.txt --output ".${value}" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 20 --beta 14 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 16 --fdr 17 --sig_threshold 0.05 --effect_allele 0 --tissue "Naive T cells versus T helper cells" --type enhancer-promoter --cis_trans cis --evidence_weight 2
done
done


for a in $datasets/promoter-enhancer/Mumbach2017/ng.3963-S5.2.tsv 
do
my_file=$(basename $a)
for value in x y
do
proc=${my_file%.tsv}_${value}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#sed -i -r 's/ \(\+[0-9]+\)//g' $processed_giggle_output 
#sed -i -r 's/ \(\-[0-9]+\)//g' $processed_giggle_output 

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table.py --tab ${proc}.lookup --input_path $a \
--tab_id Resource_index.txt --output ".${value}" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 20 --beta 14 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 16 --fdr 17 --sig_threshold 0.05 --effect_allele 0 --tissue "Naive T cells versus T regulatory cells" --type enhancer-promoter --cis_trans cis --evidence_weight 2
done
done

for a in $datasets/promoter-enhancer/Mumbach2017/ng.3963-S5.3.tsv 
do
my_file=$(basename $a)
for value in x y
do
proc=${my_file%.tsv}_${value}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#sed -i -r 's/ \(\+[0-9]+\)//g' $processed_giggle_output 
#sed -i -r 's/ \(\-[0-9]+\)//g' $processed_giggle_output 

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --input_path $a \
--tab_id Resource_index.txt --output ".${value}" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 21 --beta 14 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 16 --fdr 17 --sig_threshold 0.05 --effect_allele 0 --tissue "T helper cells versus T regulatory cells" --type enhancer-promoter --cis_trans cis --evidence_weight 2
done
done

for a in $results/coloc/transeqtl_top_hits_eqtlgen.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 9 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id_string "eQTLgen_eQTL" --study_id_column "eQTLgen" --output "_trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 9 --index_locus 1 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 2 --fdr 15 --sig_threshold 0.05 --effect_allele 7 --tissue "whole blood" --sample_size "NrSamples" --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Alasoo2018/41588_2018_46_MOESM3_ESM.5.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 14 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 14 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 14 --gene 3 --index_locus 21 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 16 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.90 --effect_allele 0 --tissue "macrophages" --sample_size 85 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Andiappan2015/ncomms8971-s2.1.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 6 --index_locus 13 --index_locus_string "" --beta 1 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "neutrophils" --sample_size 114 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.1.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.2.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "sQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.3.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "aseQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.4.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.5.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "sQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Battle2013/Supplemental_Data_1.6.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 1 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 922 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Bonder2017/exon_level_eQTLs_independent_effects.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 11 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp


#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 11 --index_locus 14 --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue "whole blood" --sample_size 3841 --type "exon eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Bonder2017/exon-ratio_level_eQTLs_independent_effects.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 11 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 11 --index_locus 14 --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue "whole blood" --sample_size 3841 --type "exon-ratio eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Bonder2017/gene_level_eQTLs.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 5 --index_locus 23 --index_locus_string "" --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 22 --sig_threshold 0.05 --effect_allele 10 --tissue "whole blood" --sample_size 3841 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Bonder2017/gene_level_eQTLs_independent_effects_interactions.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 5 --index_locus 14 --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue "whole blood" --sample_size 3841 --type "context specific eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Bonder2017/polyA-ratio_level_eQTLs_independent_effects.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 13 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 13 --index_locus 14 --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 12 --sig_threshold 0.05 --effect_allele 9 --tissue "whole blood" --sample_size 3841 --type "polyA ratio level eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Bonder2017/ng.3721-S7.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 0 --index_locus 6 --index_locus_string "" --beta 4 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 3 --tissue "whole blood" --sample_size 3841 --type "mQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

done

for a in $datasets/eQTL/Bonder2017/ng.3721-S9.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 6 --index_locus_string "" --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "whole blood" --sample_size 3841 --type "eQTM" --cis_trans trans --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Brown2017
do
for var in fat blood skin lcl
do
my_file=ng.3979-S3_${var}.txt

#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $''  \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 2 --index_locus 7 --index_locus_string "" --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 436 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/eQTL/Chen2016
do
for var in mono neut tcel
do
#Filter Allele-specific results to only significant ones
my_file=${var}_gene_WASP_ASE_all.txt	
#cat $a/$my_file | awk -v OFS="\t" '$7 < 0.05' >${my_file%.txt}.sig

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/' ${my_file%.txt}.sig  >${my_file%.txt}.sig.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.sig.short --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $'_' --index_e 2 \
--head 0 --rsid 2 --gene 3 --index_locus 9 --index_locus_string "" --beta 5 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 7 --sig_threshold 0.05 --effect_allele 1 --tissue $var --sample_size 197 --type "aseQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/eQTL/Dimas2009/GENCORD2_EQTL_FDR10_F_183
do
my_file=$(basename $a)
#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9_]+/\1/' $a >${my_file}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "fibroblast" --sample_size 197 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Dimas2009/GENCORD2_EQTL_FDR10_L_185
do
my_file=$(basename $a)
#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9_]+/\1/' $a >${my_file}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 197 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Dimas2009/GENCORD2_EQTL_FDR10_T_186
do
my_file=$(basename $a)
#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9_]+/\1/' $a >${my_file}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "T cells" --sample_size 197 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Dimas2009/GENCORD2_MQTL_FDR10_F_107
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 0 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "fibroblast" --sample_size 197 --type "mQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Dimas2009/GENCORD2_MQTL_FDR10_L_111
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 0 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 197 --type "mQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Dimas2009/GENCORD2_MQTL_FDR10_T_66
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 0 --index_locus 11 --index_locus_string "" --beta 8 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "T cells" --sample_size 197 --type "mQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Ding2010/NN57subj_p1E-5_annotate_allcis_1Mb.tbl
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 17 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 17 --index_locus 23 --index_locus_string "" --beta 7 --beta_se 8 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "healthy skin" --sample_size 55 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Ding2010/PN53subj_p1E-5_annotate_allcis_1Mb.tbl
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 17 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 17 --index_locus 23 --index_locus_string "" --beta 7 --beta_se 8 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "non-lesional skin in psorasis" --sample_size 55 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Ding2010/PP53subj_p1E-5_annotate_allcis_1Mb.tbl
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 2 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 17 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 17 --index_locus 23 --index_locus_string "" --beta 7 --beta_se 8 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "lesional skin in psorasis" --sample_size 55 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Fairfax2012/ng.2205-S2.2.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 6 --index_locus 27 --index_locus_string "" --beta 14 --beta_se 15 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 13 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "monocyte" --sample_size 288 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Fairfax2012/ng.2205-S2.3.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 6 --index_locus 27 --index_locus_string "" --beta 14 --beta_se 15 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 13 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "B cells" --sample_size 288 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Fairfax2012/ng.2205-S3.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 6 --index_locus 18 --index_locus_string "" --beta 14 --beta_se 15 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 13 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "monocyte" --sample_size 288 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Fairfax2012/ng.2205-S3.2.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 6 --index_locus 18 --index_locus_string "" --beta 14 --beta_se 15 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 13 --fdr 0 --sig_threshold 0.05 --effect_allele 3 --tissue "B cells" --sample_size 288 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Ishigaki2017/eQTL_exon_level/permutation
do
for var in B CD4 CD8 Mono NK PB
do
my_file=${var}_permutation_q_value_0.5.txt	

#Remove dots after ENSG IDs.
sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a/$my_file >${my_file%.txt}_exon.short 

python $scripts/sync_ids.py --tab ${my_file%.txt}_exon.short --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}_exon.temp

python $scripts/do_lookups.py --tab ${my_file}_exon.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 3 --delim $'\t' --delim_c $'' >${my_file}_exon.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}_exon.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 5 --index_locus 12 --index_locus_string "" --beta 9 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 105 --type "exon eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/eQTL/Ishigaki2017/eQTL_gene_level/permutation
do
for var in B CD4 CD8 Mono NK PB
do
my_file=${var}_permutation_q_value_0.5.txt	

#Remove dots after ENSG IDs.
sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a/$my_file >${my_file%.txt}_gene.short 

python $scripts/sync_ids.py --tab ${my_file%.txt}_gene.short --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}_gene.temp

python $scripts/do_lookups.py --tab ${my_file}_gene.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 3 --delim $'\t' --delim_c $'' >${my_file}_gene.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}_gene.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 4 --index_locus 11 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 105 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/eQTL/Kasela2017/journal.pgen.1006643.s010.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 2 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 2 --my_id 13 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 2 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 2 --rsid 2 --gene 13 --index_locus 18 --index_locus_string "" --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 15 --sig_threshold 0.05 --effect_allele 10 --tissue "CD4+ T cells" --sample_size 313 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Kasela2017/journal.pgen.1006643.s010.2.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 2 --my_id 2 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 2 --my_id 13 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 2 \
#--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 2 --rsid 2 --gene 13 --index_locus 18 --index_locus_string "" --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 1 --fdr 15 --sig_threshold 0.05 --effect_allele 10 --tissue "CD8+ T cells" --sample_size 313 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Kim2014/ncomms6236-s4.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 4 --gene 1 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 9 --sig_threshold 0.05 --effect_allele 0 --tissue "monocytes" --sample_size 185 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Kim2014/ncomms6236-s5.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 4 --gene 1 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 9 --sig_threshold 0.05 --effect_allele 0 --tissue "monocytes+LPS" --sample_size 185 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Lappalainen2013/EUR373.exon.cis.FDR5.all.rs137.txt 
do
my_file=$(basename $a)

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a >${my_file%.txt}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 13 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 462 --type "exon eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/eQTL/Lappalainen2013/EUR373.gene.cis.FDR5.all.rs137.txt 
do
my_file=$(basename $a)

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a >${my_file%.txt}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 13 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 462 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Lappalainen2013/EUR373.trratio.cis.FDR5.all.rs137.txt
do
my_file=$(basename $a)

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a >${my_file%.txt}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 13 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 462 --type "transcript ratio eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/LloydJones2017/cage_association_summary_statistics.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 16 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 16 --index_locus 17 --index_locus_string "" --beta 8 --beta_se 9 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole blood" --sample_size 2765 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Murphy2010/ddq392supp_table1.1.tsv
do
my_file=$(basename $a)
#Obtain only gene names in the first column
#awk -F "\t" -v OFS="\t" '{res=$1; split (res, a, " "); print a[1], $2, $3, $4}' $a >${my_file}.filtered

#python $scripts/sync_ids.py --tab ${my_file}.filtered --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 1 --index_locus 5 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "CD4+ T cells" --sample_size 200 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Naranbhai2015/ncomms8545-s2.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 8 --index_locus 9 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 5 --sig_threshold 0.05 --effect_allele 0 --tissue "neutrophils" --sample_size 101 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Naranbhai2015/ncomms8545-s4.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 7 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 7 --index_locus 8 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 5 --sig_threshold 0.05 --effect_allele 0 --tissue "neutrophils" --sample_size 101 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Nedelec2016/1-s2.0-S0092867416313071-mmc5.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 8 --index_locus 26 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 4 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "macrophage" --sample_size 91 --type "eQTL with positive selection" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Pala2017/eQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k
do
my_file=$(basename $a)

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a >${my_file%.txt}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 14 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 14 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 14 --gene 6 --index_locus 15 --index_locus_string "" --beta 9 --beta_se 10 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole peripheral blood" --sample_size 606 --type " eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Pala2017/isoQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k
do
my_file=$(basename $a)

#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9ENST_\.]+/\1/g' $a >${my_file%.txt}.short 

#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 14 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 14 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 14 --gene 6 --index_locus 15 --index_locus_string "" --beta 9 --beta_se 10 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole peripheral blood" --sample_size 606 --type " isoQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Raj2014/tableS11_meta_monocytes_cis_fdr05.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 8 --index_locus_string "" --beta 5 --beta_se 6 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "monocytes" --sample_size 401 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Raj2014/tableS12_meta_cd4T_cis_fdr05.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 8 --index_locus_string "" --beta 5 --beta_se 6 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "CD4+ T cells" --sample_size 407 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Quach2016/1-s2.0-S009286741631306X-mmc2.1.tsv
do
my_file=$(basename $a)

#sed 's/Monocyte_\t/Monocyte\t/' $a >${my_file}.filtered

#python $scripts/sync_ids.py --tab ${my_file}.filtered --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 4 --index_locus 22 --index_locus_string "" --beta 14 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 12 --fdr 13 --sig_threshold 0.05 --effect_allele 0 --tissue "Condition_desc" --sample_size 200 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Quach2016/1-s2.0-S009286741631306X-mmc4.3.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 7 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 4 --gene 7 --index_locus 25 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 16 \
--pvalue 18 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "monocyte" --sample_size 200 --type "eQTL with positive selection" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Walsh2016/13059_2016_948_MOESM1_ESM.txt
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 5 --index_locus 10 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 7 --fdr 0 --sig_threshold 0.05 --effect_allele 8 --tissue "whole blood" --sample_size 377 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Walsh2016/13059_2016_948_MOESM3_ESM.txt
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
--head 1 --my_id 6 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 6 --index_locus 7 --index_locus_string "" --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 5 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 377 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Xia2012/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt
do
my_file=$(basename $a)
#Prefilter results to 0.05 Q-value cutoff.
#(head -1 $a; awk -v OFS="\t" '$6 <0.05 {print $0}' $a) >${a%.txt}.filtered

#python $scripts/sync_ids.py --tab ${a%.txt}.filtered --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 5 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 6 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 1397 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Xia2012/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt
do
my_file=$(basename $a)

#(head -1 $a; awk -v OFS="\t" '$6 <0.05 {print $0}' $a) >${a%.txt}.filtered

#python $scripts/sync_ids.py --tab ${a%.txt}.filtered --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 5 --index_locus 7 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 6 --sig_threshold 0.05 --effect_allele 0 --tissue "LCL" --sample_size 1397 --type "eQTL" --cis_trans trans --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Yao2017/1-s2.0-S0002929717300708-mmc5.1.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 10 --index_locus_string "" --beta 3 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 7 \
--pvalue 4 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "whole blood" --sample_size 5257 --type "eQTL" --cis_trans trans-mediation --evidence_weight 2  --enable_counts 1
done

for a in $datasets/eQTL/Ye2014/tableS7-cis_eQTL_meta.8.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 20 --index_locus_string "" --beta 4 --beta_se 5 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 3 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "T cells" --sample_size 348 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/hQTL/Pelikan2018
do
for var in H3K4me1 H3K27ac
do
my_file=41467_2018_5328_MOESM4_ESM.1_${var}.tsv
#python $scripts/sync_ids.py --tab $a/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 4 --index_locus 13 --index_locus_string "" --beta 6 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 5 --fdr 0 --sig_threshold 0.05 --effect_allele 7 --tissue "LCL ${var}" --sample_size 462 --type "hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/hQTL/Pelikan2018/41467_2018_5328_MOESM10_ESM.1.tsv 
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 5 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 5 --gene 1 --index_locus 22 --index_locus_string "" --beta 10 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 7 --sig_threshold 0.05 --effect_allele 0 --tissue "Histone type" --sample_size 25 --type "hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/hQTL/Pelikan2018/41467_2018_5328_MOESM9_ESM.1.tsv
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 6 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 6 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 6 --gene 2 --index_locus 23 --index_locus_string "" --beta 11 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 8 --sig_threshold 0.05 --effect_allele 0 --tissue "Histone type" --sample_size 25 --type "hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/mQTL/GoDMC/snps_36cohorts16_rsid_vs_interval_r2_0.2_1k_significant_mr_genes
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 29 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string "GoDMC_mQTLs" --study_id_column "GoDMC" \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 29 --index_locus 31 --index_locus_string "" --beta 8 --beta_se 20 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 21 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole blood" --sample_size "TotalSampleSize" --type "mQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1 
done


for a in $datasets/network/Fagny2017
do
for var in skin whole_blood
do
my_file=${var}.qtls
#Remove dots after ENSG IDs.
#sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/g' $a/$my_file >${my_file%.txt}.short 
#python $scripts/sync_ids.py --tab ${my_file%.txt}.short --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 4 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 4 --index_locus 15 --index_locus_string "" --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 10 --sig_threshold 0.05 --effect_allele 0 --tissue $var --sample_size 310 --type "eQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/pQTL/Emilsson2018/aaq1327_Excel_tables.10.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 14 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 14 --index_locus 20 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 19 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue 13 --sample_size 15 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/pQTL/Emilsson2018/aaq1327_Excel_tables.10.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 14 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 14 --index_locus 20 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 20 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue 13 --sample_size 15 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/pQTL/Emilsson2018/aaq1327_Excel_tables.12.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".cis" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 2 --index_locus 10 --index_locus_string "" --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole blood" --sample_size 6 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/pQTL/Emilsson2018/aaq1327_Excel_tables.12.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 3 --index_locus 10 --index_locus_string "" --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 8 --fdr 0 --sig_threshold 0.05 --effect_allele 4 --tissue "whole blood" --sample_size 6 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done

for a in $datasets/pQTL/Emilsson2018/aaq1327_Excel_tables.9.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 1 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 5 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 1 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 1 --gene 5 --index_locus 17 --index_locus_string "" --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 9 --fdr 0 --sig_threshold 0.05 --effect_allele 2 --tissue "whole blood" --sample_size 6 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/pQTL/Suhre2017
do
my_file="all_suhre2017.txt"
#Include only significant p-values for Suhre et al.
#awk -F "\t" -v OFS="\t" '$11 < 10e-11 {print $0}' $a/$my_file >all_suhre2017.txt_vs_interval_r2_0.2_1k.sig

#python $scripts/sync_ids.py --tab $a/all_suhre2017.txt_vs_interval_r2_0.2_1k.sig --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 2 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 2 --index_locus 12 --index_locus_string "" --beta 9 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 11 --fdr 0 --sig_threshold 0.05 --effect_allele 6 --tissue "whole blood" --sample_size 8 --type "pQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done


for a in $datasets/RNAs/miRNA/Ziebarth2012/experimentally_supported_miRSNP_human.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 1 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "all" --sample_size "" --type "miRNA target sites" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/RNAs/miRNA/Ziebarth2012/Genes_associated_with_human_diseases_traits.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 5 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 5 --gene 1 --index_locus 8 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "all" --sample_size "" --type "miRNA target sites" --cis_trans "" --evidence_weight 3  --enable_counts 1
done

for a in $datasets/RNAs/miRNA/Ziebarth2012/SNPs_and_indels_in_miRNA_seeds_human.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 5 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 5 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 5 --gene 1 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "all" --sample_size "" --type "miRNA seed regions" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/RNAs/miRNA/Ziebarth2012/target_miRSNP_human_CLASH.txt
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 7 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 3 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 7 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 7 --gene 3 --index_locus 15 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "all" --sample_size "" --type "miRNA target sites" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM9_ESM.tsv
do
my_file=$(basename $a)

#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 13 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 1 --my_id 1 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 13 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id $datasets/Resource/Resource_index.txt --output ".trans" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 13 --gene 1 --index_locus 15 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 5 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "CD4+ T cells" --sample_size "20" --type "COGS" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $datasets/promoter-enhancer/Freire-Pritchett2017/elife-21926-supp1-v2.txt
do
for var in bait pir
do
my_file=elife-21926-supp1-v2.${var}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $',' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$var" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 17 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "hESC" --type $var --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 


for a in $datasets/promoter-enhancer/Javierre2016/PCHiC_vs_rCHiC_peak_matrix.tsv
do
for var in bait oe
do
my_file=PCHiC_vs_rCHiC_peak_matrix.${var}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 9 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$var" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 9 --index_locus 22 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "whole blood" --type $var --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 

for a in $datasets/promoter-enhancer/Mifsud2015/TS5_CD34_promoter-promoter_significant_interactions.txt
do
for var in x 
do
my_file=TS5_CD34_promoter-promoter_significant_interactions.${var}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 9 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".${var}" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 9 --index_locus 14 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "CD34+ hematopoietic cells" --type $var --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 

for a in $datasets/promoter-enhancer/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt
do
for var in x 
do
my_file=TS5_GM12878_promoter-promoter_significant_interactions.${var}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 9 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$var" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 9 --index_locus 14 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "GM12878" --type $var --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 

for var in CD34 GM12878
do
for type in x y
do
a=$datasets/promoter-enhancer/Mifsud2015/TS5_${var}_promoter-other_significant_interactions.txt
my_file=TS5_${var}_promoter-other_significant_interactions.${type}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 9 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup
python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$type" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 9 --index_locus 12 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --type $type --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.1.tsv
do
my_file=stem_2032_mmc5.1.a.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 9 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC" --type enhancer --cis_trans "" --evidence_weight 2  --enable_counts 1
done 

for a in $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.1.tsv
do
my_file=stem_2032_mmc5.1.b.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' \
#| python $scripts/sync_ids.py --db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper \
#--head 0 --my_id 8 --delim $'\t' --delim_c $';' --convert upper >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes" --delim $'\t' --delim_c $';' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 8 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 9 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes" --type enhancer --cis_trans "" --evidence_weight 2  --enable_counts 1
done 





