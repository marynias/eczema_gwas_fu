#!/bin/bash
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/ukbiobank
bgenix=/panfs/panasas01/sscm/qh18484/bin/bgenix/bin
ukbiobank=/panfs/panasas01/sscm/qh18484/ukbiobank_genetic/2017-07-04/data

HOME=/panfs/panasas01/sscm/qh18484
gwas=$HOME/data/gwas/paternoster2015
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel
utils=$HOME/bin/eczema_gwas_fu/utils
cd $analysis
#Standard exclusions:
ls $ukbiobank/derived/standard_exclusions/data.combined_recommended.txt
#Exclusion based on kinship:
ls $ukbiobank/derived/relateds_exclusions/data.minimal_relateds.txt
ls $ukbiobank/derived/relateds_exclusions/data.highly_relateds.txt
#Ancestry restrictions
ls $ukbiobank/derived/ancestry/data.europeans.bolt-lmm_only.txt
#include European sample
#Bgen files for Europeans with standard exclusions, and non-HRC markers removed.
#Need to filter those for kinship, in addition.
ls $ukbiobank/dosage_bgen
#It has been copied into the folder below to allwo analysis on the cluster.
ls /panfs/panasas01/sscm/qh18484/data/ukbiobank/dosage_bgen

#Remove the related individuals in both the lists minimal and highly relateds lists.
table_to_remove="ukbiobank_to_remove_relateds_all.txt" 
cat $ukbiobank/derived/relateds_exclusions/data.minimal_relateds.txt | cut -d" " -f1 >$table_to_remove
cat $ukbiobank/derived/relateds_exclusions/data.highly_relateds.txt | cut -d" " -f1 >>$table_to_remove

#Extract a big interval around our index SNPs - 3Mbp in each direction.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
start="$((pos - 1500000))"
end="$((pos + 1500000))"
printf -v chrom2 "%02d" $chrom
$bgenix/bgenix -g $ukbiobank/dosage_bgen/data.chr${chrom2}.bgen \
-incl-range ${chrom2}:$start-$end -vcf | gzip -c > chr${chrom}_${rsid}_1500000_ukbb.vcf.gz
python $utils/create_batch_scripts.py $utils/vcf_no_mono.sh sub_vcf_no_mono_chr${chrom}_${rsid}_1500000.sh \
my_file=chr${chrom}_${rsid}_1500000_ukbb.vcf.gz table=ukbiobank_to_remove_relateds_all.txt
qsub sub_vcf_no_mono_chr${chrom}_${rsid}_1500000.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt

while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
start="$((pos - 1500000))"
end="$((pos + 1500000))"
printf -v chrom2 "%02d" $chrom
$bgenix/bgenix -g $ukbiobank/dosage_bgen/data.chr${chrom2}.bgen \
-incl-range ${chrom2}:$start-$end  > chr${chrom}_${rsid}_1500000_ukbb.bgen
$bgenix/bgenix -g chr${chrom}_${rsid}_1500000_ukbb.bgen -index
done < $gwas/paternoster_2015_index_snps_sorted.txt

#The sample size is too big to generate the LD with current methods. 
#Use LD store instead. 

#Generate a VCF file with just one individual for the entire genome, 
#to extract positions on loci
for bgen in $ukbiobank/dosage_bgen/*.bgen
do
bgen_file=$(basename "$bgen")	
qctool_v2.0.1 -g $bgen -og ${bgen_file%.bgen}.vcf -s $ukbiobank/sample-stats/${bgen_file%.bgen}.sample -incl-samples-where 'ID_1 = 1925288'
done

#Conversion removes padding 0s for chromosome numbers, for comparison with Eczema GWAS. 
for chrom in {01..22}
do
cat $analysis/data.chr${chrom}.vcf | grep -v "#" | awk -v OFS="\t" '{print $1, $3, '0', $2, $4, $5}' | sed 's/^0//' >$analysis/diagnostics/data.chr${chrom}.bim
done

for chrom in {1..9}
do
cat $analysis/diagnostics/data.chr0${chrom}.bim | sed "s/\(rs[[:digit:]]\+\),$chrom:[^[:space:]]*/\1/" >$analysis/diagnostics/data.chr0${chrom}_modified.bim
done

for chrom in {10..22}
do
cat $analysis/diagnostics/data.chr${chrom}.bim | sed "s/\(rs[[:digit:]]\+\),$chrom:[^[:space:]]*/\1/" >$analysis/diagnostics/data.chr${chrom}_modified.bim
done

#Concatenate all bim files into 1.
touch $analysis/diagnostics/data.all.bim
for chrom in {01..22}
do
cat $analysis/diagnostics/data.chr${chrom}_modified.bim >>$analysis/diagnostics/data.all.bim
done

#Marker QC diagnosis between GWAS and SNP dataset (both with 23&me and without it).
cd $analysis/diagnostics
for my_bim in data*modified.bim data.all.bim
do
python $utils/compare_bim_files_pos.py $gwas/results.published.bim $my_bim >results.published_vs_${my_bim%.bim}.log
done
cd $analysis

#Generate GWAS table fixed in the following way:
#Remove variants where we have different alleles at the same positions, identified by the same rsids.
#Change the rsids of variants to match that in the ukbiobank dataset in cases where we have the same position, same chromosome, same alleles but different rsids.
#Remove variants present in GWAS file but not in the reference.
#Following this procedure we can match the new table using rsids in any other datasets imputed to UKBiobank.
python $utils/harmonize_rsid.py $analysis/diagnostics/data.all.bim $gwas/results.euro.bim >$gwas/results.euro.ukbiobank.bim
python $utils/update_rsid.py --bim $gwas/results.euro.ukbiobank.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$gwas/results.euro.pval.ukbiobank

#Harmonize the effect sizes so that they are relevant to the same reference allele as in the UKBiobank file. 
python $utils/harmonize_beta.py --tab $gwas/results.euro.pval.ukbiobank --ref $analysis/diagnostics/data.all.bim --header_tab Y --header_ref N \
--rsid_tab 12 --rsid_ref 2 --effect_tab 4 --alt_tab 5 --effect_ref 5 --beta_tab 8 --zscore_tab 10 --out $gwas/results.euro.pval.ukbiobank.harmonized_beta

#Generate GWAS table fixed in the following way:
#Keep all the variants and update the ID column with rsids of those found to be matching in the ref panel.
python $utils/update_rsid_all.py --bim  $analysis/diagnostics/data.all.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 --rsid 1 >$gwas/results.euro.pval.all.ukbiobank


#Generate a subset of data just with non-British Europeans to be used as a temporary reference panel.
#Take all the individuals above and filter for those in the file: 
cut -f1 -d" " $ukbiobank/derived/ancestry/data.white_british.txt >white_british_to_remove.txt

for a in *_no_rels_no_mono.vcf.gz
do
python $utils/create_batch_scripts.py $utils/vcf_no_mono2.sh sub_vcf_no_mono_${a%.vcf.gz}.sh \
my_file=$a table=white_british_to_remove.txt output=${a%.vcf.gz}_no_brits.vcf.gz
qsub sub_vcf_no_mono_${a%.vcf.gz}.sh
done

#Subset to random 10 thousand individuals among the European panel. 
#First, create a list of all the individuals contained in a random file used.
vcf-query -l chr2_rs10199605_1500000_ukbb_no_rels_no_mono_no_brits.vcf.gz >european_sample_names.txt
shuf -n 10000 european_sample_names.txt > european_sample_names_10k.txt

for a in *_no_mono_no_brits.vcf.gz
do
python $utils/create_batch_scripts.py $utils/vcf_no_mono_keep.sh sub_vcf_no_mono_${a%.vcf.gz}.sh \
my_file=$a table=european_sample_names_10k.txt output=${a%.vcf.gz}_10k.vcf.gz
qsub sub_vcf_no_mono_${a%.vcf.gz}.sh
done

#Prefilter to remove SNPs with more missing data rate than 5% 
for a in *_no_mono_no_brits_10k.vcf.gz
do
python $utils/create_batch_scripts.py $utils/vcf_no_mono_missing.sh sub_vcf_no_mono_${a%.vcf.gz}.sh \
my_file=$a maxm=0.95 output=${a%.vcf.gz}_5_perc_missing.vcf.gz
qsub sub_vcf_no_mono_${a%.vcf.gz}.sh
done

for a in *_no_mono_no_brits.vcf.gz
do
python $utils/create_batch_scripts.py $utils/vcf_no_mono_missing.sh sub_vcf_no_mono_${a%.vcf.gz}.sh \
my_file=$a maxm=0.95 output=${a%.vcf.gz}_5_perc_missing.vcf.gz
qsub sub_vcf_no_mono_${a%.vcf.gz}.sh
done


#Rename all the files for index SNP rs145809981 to rs41293864 see explaination in eczema_index_SNP_check.txt)
for a in *rs145809981*
do
b=$(echo $a | sed 's/rs145809981/rs41293864/g')
mv $a $b
done
