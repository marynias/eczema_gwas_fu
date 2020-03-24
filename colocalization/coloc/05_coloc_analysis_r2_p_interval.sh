#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/r2_interval
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
scripts_ref=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel
onek=$HOME/analysis/bayesian_fm/RefPanel/1kGenomes
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
ukbiobank=$HOME/analysis/bayesian_fm/RefPanel/ukbiobank
utils=$HOME/bin/eczema_gwas_fu/utils
#Repeat the colocalisation analysis on eQTLgen, TwinsUK and Sun 2018 using intervals defined according to r2 and p-values.

#Scenario 1: r2 >0.2 in relation to focal SNP and within 1Mbp interval centred on it. r2 values based on 1K EUR.
#Scenario 2: r2 >0.2 in relation to focal SNP and within 1Mbp interval centred on it. r2 values based on UKBiobank EUR.
#Scenario 3: p2 values higher than 1e-5 and and within 1Mbp interval centred on focal SNP.

cd $analysis
#Define GWAS SNPs to be included in the analysis based on the criteria above.
#Scenario 1
for my_ld in $onek/*_3Mbp_1kEUR.ld.gz
do
snp=$(basename $my_ld | cut -d"_" -f1)
qsub -v my_ld_file=$my_ld,my_r2=0.2,my_dataset="1k",my_snp=$snp $scripts_ref/sub_filter_by_r2.sh 
done
#Scenario 2
for my_ld in $ukbiobank/*_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.ld.gz 
do
snp=$(basename $my_ld | cut -d"_" -f2)
qsub -v my_ld_file=$my_ld,my_r2=0.2,my_dataset="ukbiobank",my_snp=$snp $scripts_ref/sub_filter_by_r2.sh 
done
#Scenario 3
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
my_start=$(expr $pos - 1500000)
my_end=$(expr $pos + 1500000)
cat $gwas/results.euro.pval.1k | awk "\$2 == $chrom {print \$0}" | awk "\$3 >= $my_start && \$3 <= $my_end {print \$0}" | awk -v OFS='\t' "\$11 <= 0.00005 {print \$3, \$11}" >${rsid}_${pos}_15000_p0.05.gwas
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Limit ourselves just to 3Mbp.
while read line
do
set -- $line
rsid=$1
pos=$3
my_start=$(expr $pos - 500000)
my_end=$(expr $pos + 500000)
cat ${rsid}_${pos}_15000_0.2.ukbiobank | awk -v OFS='\t' "\$1 >= $my_start && \$1 <= $my_end {print \$0}" >${rsid}_${pos}_5000_0.2.ukbiobank 
cat ${rsid}_${pos}_15000_0.2.1k | awk -v OFS='\t' "\$1 >= $my_start && \$1 <= $my_end {print \$0}" >${rsid}_${pos}_5000_0.2.1k 
cat ${rsid}_${pos}_15000_p0.05.gwas | awk -v OFS='\t' "\$1 >= $my_start && \$1 <= $my_end {print \$0}" >${rsid}_${pos}_5000_p0.05.gwas 
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Sort by position
for a in rs*
do
sort -k1 $a >temp
mv temp $a
done

#Create summary information about each interval.
for my_res in *.ukbiobank
do
my_snp=$(echo $my_res | cut -d"_" -f1)
my_chrom=$(grep $my_snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_interval=$(echo $my_res | cut -d"_" -f3)
my_start=$(head -1 $my_res | cut -f1)
my_end=$(tail -1 $my_res | cut -f1)
length_interval=$(expr $my_end - $my_start)
all_snps=$(cat $gwas/results.euro.pval.1k | awk "\$2 == $my_chrom {print \$0}" | awk "\$3 >= $my_start && \$3 <= $my_end {print \$0}" | wc -l | cut -f1)
high_r2_snps=$(wc -l $my_res | cut -d" " -f1)
echo -e $my_snp'\t'$my_interval'\t'$my_start'\t'$my_end'\t'$length_interval'\t'$all_snps'\t'$high_r2_snps >>all_ukbiobank.0.2
done

for my_res in *.1k
do
my_snp=$(echo $my_res | cut -d"_" -f1)
my_chrom=$(grep $my_snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_interval=$(echo $my_res | cut -d"_" -f3)
my_start=$(head -1 $my_res | cut -f1)
my_end=$(tail -1 $my_res | cut -f1)
length_interval=$(expr $my_end - $my_start)
all_snps=$(cat $gwas/results.euro.pval.1k | awk "\$2 == $my_chrom {print \$0}" | awk "\$3 >= $my_start && \$3 <= $my_end {print \$0}" | wc -l | cut -f1)
high_r2_snps=$(wc -l $my_res | cut -d" " -f1)
echo -e $my_snp'\t'$my_interval'\t'$my_start'\t'$my_end'\t'$length_interval'\t'$all_snps'\t'$high_r2_snps >>all_1k.0.2
done

for my_res in *.gwas
do
my_snp=$(echo $my_res | cut -d"_" -f1)
my_chrom=$(grep $my_snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_interval=$(echo $my_res | cut -d"_" -f3)
my_start=$(head -1 $my_res | cut -f1)
my_end=$(tail -1 $my_res | cut -f1)
length_interval=$(expr $my_end - $my_start)
all_snps=$(cat $gwas/results.euro.pval.1k | awk "\$2 == $my_chrom {print \$0}" | awk "\$3 >= $my_start && \$3 <= $my_end {print \$0}" | wc -l | cut -f1)
high_r2_snps=$(wc -l $my_res | cut -d" " -f1)
echo -e $my_snp'\t'$my_interval'\t'$my_start'\t'$my_end'\t'$length_interval'\t'$all_snps'\t'$high_r2_snps >>all_gwas.p0.05
done

#Compare the intervals generated with the methods above.
Rscript --vanilla $scripts/concatenate_r2_results.R

#Run coloc for all the genes within 3Mbp of our index SNP but using only SNPs in the designated interval.
#Therefore use the 3Mbp-based molecular trait results but subset to only SNPs within the interval using GWAS results.





