#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/final_integration/integration
scripts=$HOME/bin/eczema_gwas_fu/final_integration
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
datasets=$HOME/working/data/Datasets
utils=$HOME/bin/eczema_gwas_fu/utils

cd $analysis

#A script to calculate value-adjusted score for eah entry in the processed table.
for a in ../lookups/processed/*.processed
do
my_filename=$(basename $a)
echo $my_filename
python $scripts/adjust_score.py $a >${a%.processed}.processed2
done 

for a in ../gene_secondary/*.processed
do
my_filename=$(basename $a)
echo $my_filename
python $scripts/adjust_score.py $a >${a%.processed}.processed2
done 

for a in ../variant/processed/*.processed
do
my_filename=$(basename $a)
echo $my_filename
python $scripts/adjust_score.py $a >${a%.processed}.processed2
done 

cp -r ../variant/processed/*.processed2 ./
cp -r ../gene_secondary/*.processed2 ./
cp -r ../lookups/processed/*.processed2 ./

Rscript --vanilla $scripts/lookups_integration.R

Rscript --vanilla $scripts/model_comparison.R

Rscript --vanilla $scripts/stacked_bar_study_types.R

Rscript --vanilla $scripts/karyoploter.R

Rscript --vanilla $scripts/merge_tables_paternoster2015
#Extract the chromosomal positions of SNPs of interest.
python $utils/filter_file_by_column.py --tab Model09_snp_ranked.txt --ref $HOME/working/data/dbSNP/SNP-Base_files_merged/dbsnptable.AFchecked.allconcat.txt.gz --header_tab Y --header_ref Y --rsid_tab 2 --rsid_ref 1 >Model09_snp_ranked.dbsnp

python $utils/filter_file_by_column.py --tab Model09_snp_ranked.txt --ref $gwas/results.euro.pval.1k --header_tab Y --header_ref Y --rsid_tab 2 --rsid_ref 12 >Model09_snp_ranked.1k
#Divide evidence for SNP and gene into appropriate folders for each locus.
mkdir snp_evidence
mkdir gene_evidence

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
mkdir -p snp_evidence/$snp
mkdir -p gene_evidence/$snp
done

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
cp ${snp}_*.snps snp_evidence/${snp}/
done

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
cp ${snp}_*.genes gene_evidence/${snp}/
done

#LD for plotting SNPs with gassocplot.
for snp in rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 
do
python $scripts/obtain_r2_matrix.py ${snp}_snps.unadjusted_Model09.interval.gassocplot \
$HOME/working/data/results/RefPanel/1kGenomes/${snp}_3Mbp_1kEUR.ld.gz ${snp}_gassocplot.ld
done

#Summary table for data resources.
cat all_combined.txt | grep -v "GWAS Catalog" | grep -v "Roadmap" | awk -v OFS="\t" '{print $1, $2, $15, $16, $17, $19}' | sort | uniq >all_combined_resources.txt

#Integration with SMR results.
#Substitute the ENSG IDs to HUGO gene names, when possible.

for a in /panfs/panasas01/sscm/qh18484/working/data/Datasets/SMR_Richardson/*
do
#Remove dots after ENSG IDs.
sed -E 's/(ENSG[0-9]+)\.[0-9]+/\1/' $a >temp
python $data_manipulation_scripts/snp_json_substitute.py --tab temp \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper --head Y  --rsid 1 --delim A >${a}.hugo
done

cd $HOME/working/data/Datasets/final_integration/SMR
Rscript --vanilla $scripts/SMR_filtering.R

#Plot figures for the paper.
Rscript --vanilla $scripts/plot_interval_stats_nonlog.R
Rscript --vanilla $scripts/new_scores_histogram.R
Rscript --vanilla $scripts/loci_comparison.R
Rscript --vanilla $scripts/heatmap_evidence_top3.R