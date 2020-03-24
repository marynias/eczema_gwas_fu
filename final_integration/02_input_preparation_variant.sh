#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/final_integration/variant
scripts=$HOME/bin/eczema_gwas_fu/final_integration
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
datasets=$HOME/working/data/Datasets
giggle_analysis=$HOME/analysis/annotation/giggle
annotation_dataset=$HOME/working/data/results/euro/annotation
utils=$HOME/bin/eczema_gwas_fu/utils

cd $analysis

#Sync 1k GWAS results.
python $scripts/sync_ids.py --tab $gwas/results.euro.pval.1k --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 12 --delim $'\t' --delim_c $'' >$gwas/results.euro.pval.1k.dbsnp

#Prepare input table based on Finemapping results.
temp=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/integration
for my_snps in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for my_soft in paintor jam finemap
do
my_file=$(basename ${my_soft}_${my_snps}.input)
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
#python $scripts/sync_ids.py --tab $temp/chr$chrom.${my_snps}_clustering/$my_file --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 2 --delim $'\t' --delim_c $'' >${my_file}.temp
python $scripts/make_input_table_ver2.py --tab ${my_file}.temp --tab_id "" --input_path "" \
--tab_id_string 1 --study_id_column $my_soft --output ${my_soft}_${my_snps} --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 2 --gene 0 --index_locus 0 --index_locus_string $my_snps --beta 0 --beta_se 0 --posterior_prob 3 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0.05 --effect_allele 0 --tissue "" --type "fine-mapping" --cis_trans "" --evidence_weight 1 --enable_counts 0
done
done

#Prepare KGGSeq pathogenicity score prediction as input.
kggseq=$HOME/analysis/annotation/KGGSeq
my_file=r2_0.2_1k.kggseq_out
python $scripts/sync_ids.py --tab $kggseq/$my_file --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 17 --delim $'\t' --delim_c $'' >${my_file}.temp
#Annotate each line with all the matching index SNPs.
python $scripts/do_lookups.py --tab ${my_file}.temp --queries $data_manipulation/interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 17 --delim $'\t' --delim_c $'' >${my_file}.lookup
#Make input table.
python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id "" --input_path "" \
--tab_id_string PRVCS  --study_id_column PRVCS --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 17 --gene 0 --index_locus 22 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 19 --bayes_factor 18 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "" --type "regulatory variant prediction" --cis_trans "" --evidence_weight 3 --enable_counts 0


for a in $datasets/promoter-enhancer/Javierre2016/ActivePromoterEnhancerLinks.tsv
do
my_file=$(basename $a)
for value in oe bait
do
proc=${my_file%.tsv}.${value}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --tab_id Resource_index.txt --input_path $a \
--study_id_column 8 --output ".${value}" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "whole blood" --type "promoter enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Wang2017/csre.tab
do
my_file=$(basename $a)
proc=${my_file%.tab}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --tab_id Resource_index.txt --input_path $a \
--study_id_column 8 --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 19 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 9 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "regulatory variant predictions" --cis_trans "" --evidence_weight 3  --enable_counts 1
done

for a in $datasets/RNAs/piRNA/Wang2018/piR_hg19_sort.bed
do
my_file=$(basename $a)
proc=${my_file%.bed}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --tab_id Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 12 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "piRNA" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/promoter-enhancer/Ziebarth2013/CTCFBSDB_all_exp_sites_Sept12_2012_human_hg19.txt
do
my_file=$(basename $a)
proc=${my_file%.txt}.bed 
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${proc}.temp

#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup

python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --tab_id Resource_index.txt --input_path $a \
--study_id_column 8 --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "CTCF sites" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/caQTL/Lander2017/*bed
do
my_file=$(basename $a)
proc=$my_file
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed
tissue=$(echo $my_file | cut -d"_" -f2)
#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${proc}.temp
#python $scripts/do_lookups.py --tab ${proc}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${proc}.lookup
echo $my_file
python $scripts/make_input_table_ver2.py --tab ${proc}.lookup --tab_id Resource_index.txt --input_path $a \
--study_id_column "" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 12 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $tissue --type "FAIRE" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $datasets/splicing/Xiong2014/hg19_spidex_sig_rsid.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 8 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 7 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "splicing" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $annotation_dataset/GWAS4D/final_result_1k_r2_0.2.txt
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 3 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 3 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id_string "GWAS4D"  --study_id_column "GWAS4D" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 3 --gene 0 --index_locus 15 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 6 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "regulatory variant prediction" --cis_trans "" --evidence_weight 3  --enable_counts 0
done


for a in $annotation_dataset/fitCons/interval_r2_0.2_1k_individual.fitCons
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 7 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 7 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id_string "fitCons"  --study_id_column "fitCons" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 7 --gene 0 --index_locus 8 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 6 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "regulatory variant prediction" --cis_trans "" --evidence_weight 3  --enable_counts 0
done

#Annotate FATHMM results with rsids:
onek=$HOME/analysis/bayesian_fm/RefPanel/1kGenomes
for a in $annotation_dataset/FATHMM-FX/interval_r2_0.2_1k_alt.tab.txt $annotation_dataset/FATHMM-FX/interval_r2_0.2_1k_ref.tab.txt 
do
python $utils/update_rsid.py --bim $onek/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim --tab $a --head Y --chrom 1 --pos 2 --ref 3 --alt 4 >${a%.txt}.rsid
done

#Merge the two together and filter out entries with no score.
cat $annotation_dataset/FATHMM-FX/*.rsid >$annotation_dataset/FATHMM-FX/interval_r2_0.2_1k.merged.rsid
cat $annotation_dataset/FATHMM-FX/interval_r2_0.2_1k.merged.rsid | awk -v OFS="\t" '$6!="--"' >temp
mv temp $annotation_dataset/FATHMM-FX/interval_r2_0.2_1k.merged.rsid

for a in $annotation_dataset/FATHMM-FX/interval_r2_0.2_1k.merged.rsid
do
my_file=$(basename $a)
#python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 8 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 8 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a \
--tab_id_string "FATHMM-FX"  --study_id_column "FATHMM-FX" --output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 8 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 6 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "regulatory variant prediction" --cis_trans "" --evidence_weight 3  --enable_counts 0
done

#Filter out significant hits showing allelic imbalance in Maurano results.
for a in $datasets/caQTL/Maurano2015/ng.3432-S5.txt
do
cat $a | awk '$12 < 0.05' >ng.3432-S5_filtered.txt
my_file=ng.3432-S5_filtered.txt

python $scripts/sync_ids.py --tab $my_file --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 14 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 12 --sig_threshold 0.05 --effect_allele 0 --tissue "all" --type "caQTL" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/caQTL/Maurano2015/ng.3432-S7.txt
do
my_file=$(basename $a)

python $scripts/sync_ids.py --tab $a --db $data_manipulation/rsid_synonyms.txt \
--head 1 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 4 --gene 0 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 5 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "caQTL" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

#GWAS catalog look-up
for a in $datasets/GWAS/gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv
do
my_file=$(basename $a)
cat $a | awk -F "\t" -v OFS="\t" '$21 !~ /;/ {print $0}' >${my_file}.temp 
cat ${my_file}.temp | awk -F "\t" -v OFS="\t" '{print $2, $3, $5, $7, $8, $21, $28, $31, $35}' >${my_file}.intermediate
awk -F "\t" -v OFS="\t" '{res=$6; split (res, a, "-"); print $1, $2"_"$4"_"$3,  $5, a[1], a[2], $7, $8, $9}' ${my_file}.intermediate >${my_file}.filtered
done


for a in gwas_catalog_v1.0.2-associations_e93_r2019-01-11.tsv.filtered
do
my_file=$a
#Convert from UTF-8 to ASCII
#cat $my_file | iconv -f utf-8 -t ascii//TRANSLIT >${my_file}.translit
#python $scripts/sync_ids.py --tab ${my_file}.translit --db $data_manipulation/rsid_synonyms.txt \
#--head 1 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 1 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

#sed -i 's/?/NA/' ${my_file}.lookup

#Eliminate datasets already analized
cat ${my_file}.lookup | grep -v "Emilsson V_Co-regulatory networks of human serum proteins link genetics to disease." | grep -v "Sun BB_Genomic atlas of the human plasma proteome._Nature" | grep -v "Suhre K_Connecting genetic risk to disease end points through the human blood plasma proteome._Nat Commun" | grep -v "atopic eczema" | grep -v "Paternoster" >${my_file}.lookup2

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup2 \
--tab_id_string 2  --study_id_column 1 --output "GWAS_catalog" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 1 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 7 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 6 --fdr 0 --sig_threshold 0 --effect_allele 5 --tissue "ALLELE" --type "GWAS Catalog" --cis_trans "" --evidence_weight 2  --enable_counts 0
done


for a in $datasets/hQTL/Chen2016
do
for var in mono neut tcel
do
for acetylation in K4ME1 K27AC
do
my_file=${var}_${acetylation}_log2rpm_peer_10_all_summary.txt
#Filter results to only significant ones
cat $a/$my_file | awk -v OFS="\t" '$7 < 0.05 {print $0}' >$a/${my_file%.txt}.sig 

python $scripts/sync_ids.py --tab $a/${my_file%.txt}.sig --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 2 --delim $' ' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 2 --delim $' ' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $'_' --index_e 2 \
--head 0 --rsid 2 --gene 0 --index_locus 10 --index_locus_string "" --beta 5 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 7 --sig_threshold 0.05 --effect_allele 0 --tissue ${var}_${acetylation} --sample_size 197 --type "hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done
done

for a in $datasets/hQTL/Chen2016
do
for var in mono neut tcel
do
for acetylation in H3K27ac H3K4me1
do
my_file=${var}.${acetylation}_peak_WASP_ASE_all.txt
#Filter results to only significant ones
cat $a/$my_file | awk -v OFS="\t" '$7 < 0.05 {print $0}' >$a/${my_file%.txt}.sig 

python $scripts/sync_ids.py --tab $a/${my_file%.txt}.sig --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 2 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 2 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --input_path $a/$my_file \
--tab_id $datasets/Resource/Resource_index.txt --output "" --delim $'\t' --delim_c $',' --delim_e $'_' --index_e 2 \
--head 0 --rsid 2 --gene 0 --index_locus 9  --index_locus_string "" --beta 5 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 4 --fdr 7 --sig_threshold 0.05 --effect_allele 0 --tissue ${var}_${acetylation} --sample_size 197 --type "ase hQTL" --cis_trans cis --evidence_weight 2  --enable_counts 1
done
done
done

for a in $giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_roadmap_sort_index_processed.bed
do
my_file=$(basename $a)
#Eliminate non-active regions
#cat $a | grep -v 15_Quies | grep -v 14_ReprPCWk | grep -v 13_ReprPC | grep -v 12_EnhBiv | grep -v 11_BivFlnk | grep -v 10_TssBiv | grep -v 9_Het >${a%.bed}.filtered.bed
#python $scripts/sync_ids.py --tab ${a%.bed}.filtered.bed --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp
#Add additional column to help distinguish signal on the same SNP from different sources
#cat ${my_file}.temp | sed -E 's:/panfs/panasas01/sscm/qh18484/bin/giggle/roadmap_sort/::g' >temp
#cat temp | awk -v OFS="\t" '$10=$8"_"$9 {print $0}' >${my_file}.temp
python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup
python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id_string 8  --study_id_column Roadmap --output "Roadmap" \
--delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 --alt_field 10 \
--head 0 --rsid 4 --gene 0 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue 9 --type "epigenome" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/TADs/Rao2014/GM12878_CTCF_orientation.bed
do
my_file=$(basename $a)
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 15 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "LCL" --type "CTCF sites" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/TADs/Harmston2017/41467_2017_524_MOESM2_ESM.bed
do
my_file=$(basename $a)

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 15 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "Converved genomic regulatory blocks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $datasets/TADs/Freire-Pritchett2017/hESC_TADs_delta2.0.bed
do
my_file=$(basename $a)

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

head $processed_giggle_output

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "hESC" --type "TAD" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $datasets/TADs/ENCODE
do
for var in RPMI7951 SKMEL5
do
my_file=ENCODE3-${var}-HindIII__hg19__genome__C-40000-iced.tads.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a/$my_file \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 11 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --type "TAD" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 

for a in $datasets/TADs/Dixon2012
do
for var in hESC IMR90
do
my_file=${var}_total.combined.domain

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}.bed_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}.bed_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a/$my_file \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue $var --type "TAD" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 

for a in $datasets/RNAs/lncRNA/Xie2014/NONCODEv5_hg19.lncAndGene.bed
do
my_file=$(basename $a)

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 18 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "lncRNA" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/RNAs/lncRNA/Chen2012/human_hg37.bed
do
my_file=$(basename $a)

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "all" --type "lncRNA" --cis_trans "" --evidence_weight 2  --enable_counts 1
done


for a in $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM5_ESM
do
for var in bait oe
do
my_file=13059_2017_1285_MOESM5_ESM.${var}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$var" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "CD4+ T cells" --type $var --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 


for a in $datasets/promoter-enhancer/Mifsud2015/mifsud2015_supp_tab2_LCL_TADs.bed
do
my_file=$(basename $a)
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "LCL" --type "TAD" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/promoter-enhancer/Ziebarth2013/CTCFBSDB_all_exp_sites_Sept12_2012_human_hg19.txt
do
my_file="CTCFBSDB_all_exp_sites_Sept12_2012_human_hg19.bed"
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 10 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue 8 --type "CTCF" --cis_trans "" --evidence_weight 2  --enable_counts 1
done

for a in $datasets/promoter-enhancer/Wang2018
do
for my_file in 41467_2018_7746_MOESM4_ESM.txt 41467_2018_7746_MOESM6_ESM.txt
do
proc=${my_file%.txt}.bed

raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${proc}_index_processed.bed

#head $processed_giggle_output

#python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
#--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

#python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
#--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a/$my_file \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "LCL" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done 


for a in $datasets/promoter-enhancer/Wang2015/pnas.1507253112.sd01_hg19.bed
do
my_file=$(basename $a)
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output "" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "CD4+ T cell" --type "insulator" --cis_trans "" --evidence_weight 3  --enable_counts 1
done


for a in $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.1.tsv \
$datasets/promoter-enhancer/Rubin2017/ng.3935-S5.2.tsv \
$datasets/promoter-enhancer/Rubin2017/ng.3935-S5.3.tsv \
$datasets/promoter-enhancer/Rubin2017/ng.3935-S5.4.tsv 
do
my_file=$(basename $a | sed 's/tsv//')
for var in x y
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}${var}.bed_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}${var}.bed_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".$var" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "keratinocytes" --type "promoter enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc3.2_EpiSC_DNMT3A_peaks.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_DNMT3A" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC" --type "DNMT3A_peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done


for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc3.2_EpiSC_DNMT3B_peaks.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_DNMT3B" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC" --type "DNMT3B_peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc3.2_Diff_DNMT3A_peaks.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Diff_DNMT3A" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "keratinocytes" --type "DNMT3A_peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done


for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc3.2_Diff_DNMT3B_peaks.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Diff_DNMT3B" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "keratinocytes" --type "DNMT3B_peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K4me1.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_H3K4me1" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC_H3K4me1" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K4m3.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_H3K4m3" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC_H3K4m3" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done


for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K27ac.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_H3K27ac" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC_H3K27ac" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done


for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_Diff_H3K4m1.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Diff_H3K4m1" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_H3K4m1" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_Diff_H3K4m3.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Diff_H3K4m3" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_H3K4m3" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.1_Diff_H3K27ac.bed 
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Diff_H3K27ac" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_H3K27ac" --type "ChIP-seq peaks" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.a.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".EpiSC_consensus" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "EpiSC_consensus" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.b.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes_consensus" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_consensus" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.c.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes_lost" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_lost" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.d.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes_active" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_active" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done

for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.e.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes_gained" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_gained" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done


for a in $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv
do
for my_file in 1-s2.0-S1934590916301953-mmc4.2.a.bed
do
raw_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index.gz
processed_giggle_output=$giggle_analysis/interval_r2_0.2_1k_individual.bed.gz_vs_${my_file}_index_processed.bed

python $scripts/sync_ids.py --tab $processed_giggle_output --db $data_manipulation/rsid_synonyms.txt \
--head 0 --my_id 4 --delim $'\t' --delim_c $'' >${my_file}.temp

python $scripts/do_lookups.py --tab ${my_file}.temp --queries interval_r2_0.2_1k_nodups_sync.map --head 0 \
--my_id 4 --delim $'\t' --delim_c $'' >${my_file}.lookup

python $scripts/make_input_table_ver2.py --tab ${my_file}.lookup --tab_id $datasets/Resource/Resource_index.txt --input_path $a \
--output ".Keratinocytes_gained" --delim $'\t' --delim_c $',' --delim_e $',' --index_e 0 \
--head 0 --rsid 4 --gene 0 --index_locus 9 --index_locus_string "" --beta 0 --beta_se 0 --posterior_prob 0 --bayes_factor 0 --score 0 \
--pvalue 0 --fdr 0 --sig_threshold 0 --effect_allele 0 --tissue "Keratinocytes_gained" --type "enhancer" --cis_trans "" --evidence_weight 2  --enable_counts 1
done
done