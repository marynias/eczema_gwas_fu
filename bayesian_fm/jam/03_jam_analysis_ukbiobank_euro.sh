#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
JAM=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam
JAM_ANALYSIS=$HOME/analysis/bayesian_fm/jam/ukbiobank
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/jam
scripts_ref=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
utils=$HOME/bin/eczema_gwas_fu/utils
ukbb=$HOME/analysis/bayesian_fm/RefPanel/ukbiobank

cd $JAM_ANALYSIS

#Generate range of SNPs to be analyzed in each interval.
#Warning: this will output all the SNPs in a given interval, regardless of the chromosome number.
#SNPs matching the positions in the input VCF file with genotypes, which match the right chromosome for given SNP will be later used to obtain correct results.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000
do
qsub -v reference=$gwas/results.euro.pval.ukbiobank.harmonized_beta,snp=$my_snps,interval=$my_int $utils/sub_extract_snp_interval.sh
done
done

#Generate JAM input
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_input=$ukbb/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
qsub -v locus=${my_snps}.${my_int}.locus,my_vcf=$my_input,chr=${chrom},snp=$my_snps,interval=$my_int $scripts/sub_generate_jam_input.sh
done
done

#It looks like needs to prune SNPs in high LD (>0.95) and with low MAF (<0.05)
#Convert reference panel VCFs to Plink format.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_vcf=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
python $utils/create_batch_scripts.py $utils/vcf_to_plink.sh sub_vcf_to_plink_tr_${chrom}.sh input_file=$my_vcf optional="--recode transpose"
qsub sub_vcf_to_plink_tr_${chrom}.sh
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
#Add fake sex identity to tfam file, as required by PriorityPruner.
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_vcf=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
awk '{$5='1'; print $0}' ${my_vcf%.vcf*}.tfam >temp
mv temp ${my_vcf%.vcf*}.tfam
#Remove additional characters after RSID Id in the tped file.
sed -i "s/\(rs[[:digit:]]\+\)\S*/\1/" ${my_vcf%.vcf*}.tped
done

##Generate Priority Pruner input files & run it.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
python $scripts/generate_pp_snp_table2.py --tab $gwas/results.euro.pval.ukbiobank.harmonized_beta \
--proces chr${chrom}.$my_snps.$my_int.beta --proces_id 1 --out chr${chrom}.$my_snps.$my_int.snp.table \
--chrom 2 --pos 3 --ident 12 --ref 4 --alt 5 --pval 11 --forceSelect ${my_snps}
my_vcf=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
python $utils/create_batch_scripts.py $utils/priority_pruner.sh \
sub_priority_pruner_${chrom}_${my_snps}_$my_int.sh \
tfile=${my_vcf%.vcf*} \
snp_table=chr${chrom}.${my_snps}.${my_int}.snp.table \
r2=0.95 \
interval=$my_int \
min_maf=0.05 \
min_snp_call_rate=0.9 \
output=chr${chrom}_${my_snps}_${my_int}_0.95r2_0.95maf 
qsub sub_priority_pruner_${chrom}_${my_snps}_$my_int.sh
done
done

#Run JAM.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 
do
chr=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
qsub -v chrom=$chr,snp=$my_snps,interval=$my_int,r2="0.95",maf="0.95" $scripts/sub_run_jam.sh
done
done

#Carry out analysis using haploblocks generated by BigLD on 1K EUR data.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
qsub -v my_rsid=$rsid,my_chrom=$chrom,my_pos=$pos,my_bigld=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bigld,my_gwas=$gwas/results.euro.pval.ukbiobank.harmonized_beta $scripts_ref/sub_extract_snp_interval_bigld.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Generate JAM input for analysis involving haploblocks.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_input=$ukbb/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
qsub -v locus=${my_snps}.bigld.locus,my_vcf=$my_input,chr=${chrom},snp=$my_snps $scripts/sub_generate_jam_input_bigld.sh
done

#Priority Pruner
#Remember to Put back pruned SNPs with high r2 with the ones in the high PP credible interval, and colour them differently in the heatmap.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
python $scripts/generate_pp_snp_table2.py --tab $gwas/results.euro.pval.ukbiobank.harmonized_beta \
--proces chr${chrom}.$my_snps.bigld.beta --proces_id 1 --out chr${chrom}.$my_snps.bigld.snp.table \
--chrom 2 --pos 3 --ident 12 --ref 4 --alt 5 --pval 11 --forceSelect ${my_snps}
my_vcf=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
python $utils/create_batch_scripts.py $utils/priority_pruner.sh \
sub_priority_pruner_${chrom}_${my_snps}_bigld.sh \
tfile=${my_vcf%.vcf*} \
snp_table=chr${chrom}.${my_snps}.bigld.snp.table \
r2=0.95 \
interval="999999999999" \
min_maf=0.05 \
min_snp_call_rate=0.9 \
output=chr${chrom}_${my_snps}_bigld_0.95r2_0.95maf 
qsub sub_priority_pruner_${chrom}_${my_snps}_bigld.sh
done

#Run JAM.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chr=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
qsub -v chrom=$chr,snp=$my_snps,interval="bigld",r2="0.95",maf="0.95" $scripts/sub_run_jam.sh
done

###Same as above, but using gpart partitioning for generating haploblocks.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
qsub -v my_rsid=$rsid,my_chrom=$chrom,my_pos=$pos,my_bigld=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.gpart,my_gwas=$gwas/results.euro.pval.ukbiobank.harmonized_beta $scripts_ref/sub_extract_snp_interval_bigld.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Generate JAM input for analysis involving haploblocks.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
my_input=$ukbb/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
qsub -v locus=${my_snps}.locus,my_vcf=$my_input,chr=${chrom},snp=$my_snps $scripts/sub_generate_jam_input_gpart.sh
done

#Priority Pruner
#Remember to Put back pruned SNPs with high r2 with the ones in the high PP credible interval, and colour them differently in the heatmap.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
python $scripts/generate_pp_snp_table2.py --tab $gwas/results.euro.pval.ukbiobank.harmonized_beta \
--proces chr${chrom}.${my_snps}.gpart.beta --proces_id 1 --out chr${chrom}.$my_snps.gpart.snp.table \
--chrom 2 --pos 3 --ident 12 --ref 4 --alt 5 --pval 11 --forceSelect ${my_snps}
my_vcf=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz
python $utils/create_batch_scripts.py $utils/priority_pruner.sh \
sub_priority_pruner_${chrom}_${my_snps}_bigld.sh \
tfile=${my_vcf%.vcf*} \
snp_table=chr${chrom}.${my_snps}.gpart.snp.table \
r2=0.95 \
interval="999999999999" \
min_maf=0.05 \
min_snp_call_rate=0.9 \
output=chr${chrom}_${my_snps}_gpart_0.95r2_0.95maf 
qsub sub_priority_pruner_${chrom}_${my_snps}_bigld.sh
done

#Run JAM.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chr=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
qsub -v chrom=$chr,snp=$my_snps,interval="gpart",r2="0.95",maf="0.95" $scripts/sub_run_jam.sh
done


#Calculate LD statistics for each SNP pair using VCF tools.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
qsub -v my_file=${ukbb}/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.vcf.gz $utils/sub_vcf_ld.sh 
done

#Check that the LD between the non-selected SNPs in Priority Pruner and the chosen SNPs in the credible set.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 bigld gpart
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
ld=$ukbb/chr${chrom}_${my_snps}_1500000_ukbb_no_rels_no_mono_no_brits_10k_5_perc_missing.ld.gz
results=chr${chrom}_${my_snps}_${my_int}_0.95r2_0.95maf.results 
set=chr${chrom}.${my_snps}.${my_int}_credible_set.txt
r2=0.95
output=chr${chrom}.${my_snps}.${my_int}_credible_set_extended.txt
qsub -v my_ld=$ld,my_results=$results,my_set=$set,my_r2=$r2,my_output=$output $scripts/sub_generate_expanded_confidence_interval.sh 
done
done


#Divide the output with one SNP per folder
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 bigld gpart
do
cd $JAM_ANALYSIS
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
mkdir -p chr${chrom}.${my_snps}/${my_int}
cp -r chr*${my_snps}.${my_int}_* chr${chrom}.${my_snps}/${my_int}
done
done

#Rerun the pipeline for the correct SNP ID for the MICB locus - rs145809981 has been merged into rs41293864 on dbSNP. The position is missing from UKBiobank altogether.
#In addition the following index snps are missing from UKbiobank: rs12188917, rs12730935, rs2592555, rs6419573 (see eczema_index_SNP_check.txt doc for details of analysis)
#Remove the files pertaining to the 5 loci missing from UK Biobank. Therefore, we should end up with 39 loci with results in this analysis.
rm -rf *rs145809981*
rm -rf *rs41293864*
rm -rf *rs12730935*
rm -rf *rs12188917*
rm -rf *rs2592555*
rm -rf *rs6419573*

