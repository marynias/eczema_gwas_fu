#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/networks/dmgwas
scripts=$HOME/bin/eczema_gwas_fu/annotation/networks
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
utils=$HOME/bin/eczema_gwas_fu/utils

cd $analysis

#Create input file to be used for online VEGAS2 analysis to generate gene-based p-values from GWAS.
#Need to split into two batches (chrom 1-11, 12-22) as the resulting file is too big.
tail -n +2 $gwas/results.euro.pval.1k | awk -v OFS="\t" '$2 >=1 && $2 <=11 {print $12, $11}' >results.euro.pval.1k.vegas.part1
tail -n +2 $gwas/results.euro.pval.1k | awk -v OFS="\t" '$2 >=12 && $2 <=22 {print $12, $11}' >results.euro.pval.1k.vegas.part2

#Download VEGAS results.
wget http://vegas2.qimrberghofer.edu.au/download/results_euro_pval_1k_vegas_part1-6165B92C-1A58-11E9-B610-D5FCD8CB637A.v2out
wget http://vegas2.qimrberghofer.edu.au/download/results_euro_pval_1k_vegas_part2-76391AA6-1A58-11E9-B767-2BFDD8CB637A.v2out

#Join into one file (for NetWAS analyses)
tail -n +2 results_euro_pval_1k_vegas_part1-6165B92C-1A58-11E9-B610-D5FCD8CB637A.v2out >results_euro_pval_1k.vegas.out
tail -n +2 results_euro_pval_1k_vegas_part2-76391AA6-1A58-11E9-B767-2BFDD8CB637A.v2out >>results_euro_pval_1k.vegas.out
cat results_euro_pval_1k.vegas.out | awk '{$1 = "\x22"$1"\x22"; $2 = "\x22"$2"\x22"; $9 = "\x22"$9"\x22"; print}' | awk '{$4 = sprintf("%.0f", $4); $8 = sprintf("%.14f", $8); $10 = sprintf("%.14f", $10); print}' | sort -k1,2  >results_euro_pval_1k.vegas.out2 

#Substitute Gene names with HUGO identifiers
for a in results_euro_pval_1k_vegas_part1-6165B92C-1A58-11E9-B610-D5FCD8CB637A.v2out \
results_euro_pval_1k_vegas_part2-76391AA6-1A58-11E9-B767-2BFDD8CB637A.v2out
do
python $data_manipulation_scripts/snp_json_substitute.py --tab $a \
--db $data_manipulation/hugo_synonyms_ids2_filtered.hugo.upper --head Y  --rsid 2 --delim B>${a}.hugo
done

#Prepare dmGWAS input.
cat *.hugo | tail -n +2 | cut -f2,8 >results_euro_pval_1k_vegas.dmgwas

#Substitute quotes
sed -i 's/"//g' results_euro_pval_1k_vegas_part1-6165B92C-1A58-11E9-B610-D5FCD8CB637A.v2out
sed -i 's/"//g' results_euro_pval_1k_vegas_part2-76391AA6-1A58-11E9-B767-2BFDD8CB637A.v2out

#Download of protein-protein interaction dataset from STRING db.
wget https://stringdb-static.org/download/protein.links.full.v10.5/9606.protein.links.full.v10.5.txt.gz 

#Substitute the 9606. prefix
sed -i 's/9606\.//g' 9606.protein.links.full.v10.5.txt 

#Filter only for the more reliable connections with high score. Total score range 150-999.
#Filtering above 300 score, roughly 20$ of total input - 2652169 interactions left versus 11353057 initially.
#Doing this because weak interactions are only based on textmining.
cat 9606.protein.links.full.v10.5.txt | awk '$16 > 300' >9606.protein.links.full.v10.5.conf.txt

#Convert ENSP ids to gene names.
Rscript --vanilla $scripts/string_gene_names.R

#Generate interaction input as required by dmGWAS.
tail -n +2 9606.protein.links.full.v10.5.conf.ensembl | cut -f18,19 | sort >9606.protein.links.full.v10.5.conf.dmgwas

cd $analysis/Suarez-Farinas2011
Rscript --vanilla $scripts/expression_tables.R
cd $analysis

qsub -v my_script=$scripts/dmgwas.R $utils/sub_run_Rscript.sh