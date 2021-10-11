#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
JAM=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam
JAM_ANALYSIS=$HOME/analysis/bayesian_fm/jam/1k_published
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
utils=$HOME/bin/eczema_gwas_fu/utils
#Submit test jam run
cd $JAM_ANALYSIS
#qsub $JAM/sub_jam_test.sh
#The analysis below follows through from finemap_analysis.sh, and includes some of the files produced in it.


#It looks like needs to prune SNPs in high LD (>0.95) and with low MAF (<0.05)
#Convert reference panel VCFs to Plink format.
for chrom in 1 2 11
do
my_vcf=${analysis}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
python $utils/create_batch_scripts.py $utils/vcf_to_plink.sh sub_vcf_to_plink_tr_${chrom}.sh input_file=$my_vcf optional="--recode transpose"
qsub sub_vcf_to_plink_tr_${chrom}.sh
done

#Turns out that rs IDs do not sometimes match in the reference 1k file (while the positions and alleles, do!).
#In that case, substitute rsids from the our GWAS file into our tped file. 
#Add fake sex identity to tfam file, as required by PriorityPruner.
for chrom in 1 2 11
do
my_vcf=${analysis}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
python $scripts/substitute_rsids.py --tab $gwas/results.published.tsv \
--ped ${my_vcf%.vcf*}.tped --out ${my_vcf%.vcf*}.gwas.tped \
--ident 1 --chrom 2 --pos 3 >${my_vcf%.vcf*}.gwas.log
awk '{$5='1'; print $0}' ${my_vcf%.vcf*}.tfam >${my_vcf%.vcf*}.gwas.tfam
done

function jam_step1
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
#Check which chromosome our SNP is located on.
my_vcf=${analysis}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
#Script to produce input for JAM
python $scripts/generate_jam_input.py --tab $gwas/results.published.tsv \
--proces $analysis/chr${chrom}.$snp.$interval.processed --vcf $my_vcf \
--out chr${chrom}.$snp.$interval --pos 3 --ref 4 --alt 5 --beta 8 --se 9

#Generate input for Priority Pruner
python $scripts/generate_pp_snp_table.py --tab $gwas/results.published.tsv \
--proces $analysis/chr${chrom}.$snp.$interval.processed --out chr${chrom}.$snp.$interval.snp.table \
--chrom 2 --pos 3 --ident 1 --ref 4 --alt 5 --pval 10 --forceSelect $snp

vcf_name=$(basename "$my_vcf")
python $utils/create_batch_scripts.py $utils/priority_pruner.sh sub_priority_pruner_${chrom}_${snp}_$interval.sh \
tfile=${my_vcf%.vcf*}.gwas snp_table=chr${chrom}.$snp.$interval.snp.table r2=0.95 interval=$interval min_maf=0.05 \
min_snp_call_rate=0.9 output=${vcf_name%.vcf*}_${chrom}_${snp}_${interval}_0.95r2_0.95maf 
qsub sub_priority_pruner_${chrom}_${snp}_$interval.sh
}

jam_step1 rs2212434 1500000
jam_step1 rs2212434 250000
jam_step1 rs2212434 140000
jam_step1 rs2212434 50000

jam_step1 rs61813875 1500000
jam_step1 rs61813875 250000
jam_step1 rs61813875 140000
jam_step1 rs61813875 50000

jam_step1 rs112111458 1500000
jam_step1 rs112111458 250000
jam_step1 rs112111458 140000
jam_step1 rs112111458 50000

#Modify the JAM analysis script to subset to only SNPs that passed the pruning step.
#Going for r2 < 0.95 and MAF > 0.05 (when MAF > 0.01 or r2 < 0.99 
#get the error: "the leading minor of order 23 is not positive definite")

function jam_step2
{
snp=$1
interval=$2
r2=$3
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_vcf=${analysis}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
vcf_name=$(basename "$my_vcf")
Rscript $scripts/run_jam.R chr${chrom}.$snp.$interval.beta chr${chrom}.$snp.$interval.matrix ${vcf_name%.vcf*}_${chrom}_${snp}_${interval}_${r2}r2_0.95maf.results
}

jam_step2 rs2212434 1500000 0.95
jam_step2 rs2212434 250000 0.95
jam_step2 rs2212434 140000 0.95
jam_step2 rs2212434 50000 0.95

jam_step2 rs61813875 1500000 0.95
jam_step2 rs61813875 250000 0.95
jam_step2 rs61813875 140000 0.95
jam_step2 rs61813875 50000 0.95

jam_step2 rs112111458 1500000 0.95
jam_step2 rs112111458 250000 0.95
jam_step2 rs112111458 140000 0.95
jam_step2 rs112111458 50000 0.95

#Only works for very small interval - 50000. Otherwise, get the positive definite error.
#Check if decreasing r2 threshold, to < 0.8, will work.

function pruner
{
snp=$1
interval=$2
r2=$3
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_vcf=${analysis}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
vcf_name=$(basename "$my_vcf")
python $utils/create_batch_scripts.py $utils/priority_pruner.sh sub_priority_pruner_${chrom}_${snp}_$interval.sh \
tfile=${my_vcf%.vcf*}.gwas snp_table=chr${chrom}.$snp.$interval.snp.table r2=${r2} interval=$interval min_maf=0.05 \
min_snp_call_rate=0.9 output=${vcf_name%.vcf*}_${chrom}_${snp}_${interval}_${r2}r2_0.95maf 
qsub sub_priority_pruner_${chrom}_${snp}_$interval.sh
}
pruner rs2212434 1500000 0.8
pruner rs2212434 250000 0.8
pruner rs2212434 140000 0.8

jam_step2 rs2212434 1500000 0.8
jam_step2 rs2212434 250000 0.8
#This procedure results only in this interval working:
jam_step2 rs2212434 140000 0.8

##Now, check for EMSY and FLG loci. 
pruner rs61813875 1500000 0.8
pruner rs61813875 250000 0.8
pruner rs61813875 140000 0.8
pruner rs61813875 50000 0.8

pruner rs112111458 1500000 0.8
pruner rs112111458 250000 0.8
pruner rs112111458 140000 0.8
pruner rs112111458 50000 0.8

jam_step2 rs61813875 1500000 0.8
jam_step2 rs61813875 250000 0.8
jam_step2 rs61813875 140000 0.8
jam_step2 rs61813875 50000 0.8

jam_step2 rs112111458 1500000 0.8
jam_step2 rs112111458 250000 0.8
jam_step2 rs112111458 140000 0.8
jam_step2 rs112111458 50000 0.8

#None of these analyses work...


