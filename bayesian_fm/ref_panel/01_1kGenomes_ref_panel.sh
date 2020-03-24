#!/bin/bash
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
scripts=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/bayesian_fm/ref_panel
paintor=/panfs/panasas01/sscm/qh18484/bin/PAINTOR_V3.0/PAINTOR_Utilities
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
scripts_jam=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/bayesian_fm/jam
utils=$HOME/bin/eczema_gwas_fu/utils

cd $analysis

#Test run of the pipeline from Paintor to calculate LD Matrix: https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD

###Testing using filaggrin locus rs61813875	1q21.3 
#Download chromosome 1 VCF file and panel file
input_vcf=ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${input_vcf}.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

#Format locus file from raw data in results.published tsv.
python $scripts/tab2locus_paintor.py --tab $gwas/results.published.tsv \
--out $analysis/chr1.locus --chrom 2 --chr_out 1 \
--pos 3 --ident 1 --ref 4 --alt 5 --beta 8 --se 9 --indels N

#Note, the locus file needs to be sorted by position!

#Try the script calculating LD matrix
python $paintor/CalcLD_1KG_VCF.py \
--locus chr1.locus \
--reference ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--map integrated_call_samples_v3.20130502.ALL.panel \
--effect_allele A1 \
--alt_allele A0 \
--population EUR \
--Zhead Zscore \
--out_name chr1.ld \
--position pos

#Terminated early due to choosing the analysis with ambiguous SNPs below.

#A short test job
qsub $scripts/sub_CalcLD_1KG_VC_1K_short.sh

#Note that the script removes SNPs which were monomorphic in the reference panel and which yield NaN correlation values. To avoid confusion, will prefilter the input VCF file for monomorphic SNPs.

#Note that the script in its current form compares the SNPs in GWAS table and population reference VCF file by position and SNP id, requiring the use of the same genome reference build! For 1k and eczema GWAS this holds true - both use GRCh37.p13	but that may not be the case for other reference panels.

#Note that SNPs for which more than 2 alleles are found (multimorphic) in the ref panel are ignored in calculating the output. 

#The script outputs two files: ld with the LD matrix and processed which filters their input locus file for SNPs which were included in the LD matrix calculation.

#First filter out non-European samples from the input VCF file.

#Capture the individual sample names to be filtered in a single variable. 
table=integrated_call_samples_v3.20130502.ALL.panel
to_remove=$(awk '($3 != "EUR") {print "--remove-indv " $1}' $table | tr '\n' ' ')

#Screen the VCF file to eliminate monomorphic SNPs in the ref panel.
#Option --recode is essential for the program to output any variants!
qsub -t 1 $scripts/sub_vcf_no_mono_1k.sh

#Create a population mapping table with only European individuals
awk '{if (NR==1) print $0 ; else if ($3 == "EUR") print $0}' $table >${table%.ALL.panel}.EUR.panel

#Generate locus file for 3 Mbp around the FLG locus index SNP. 
Rscript --vanilla $scripts/extract_snp_interval.R chr1.locus rs61813875 1500000

#Run the LD matrix generation script with the subset
python $utils/create_batch_scripts.py $scripts/sub_CalcLD_1KG_VC_1K.sh  \
sub_CalcLD_rs61813875_1500000.sh my_file=chr1.rs61813875.1500000.locus chrom=1
qsub sub_CalcLD_rs61813875_1500000.sh

#Obtained correct results - some loci from GWAS will not be available due to them containing more than 2 alleles in reference panel, and being filtered out.

#Run the analysis above for chromosome 2 and index SNP rs112111458 (CD207) and chromosome 11 and index snp rs2212434 (EMSY/C11orf30). 
function get_chrom_ref 
{
chrom=$1
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
python $scripts/tab2locus_paintor.py --tab $gwas/results.published.tsv \
--out $analysis/chr${chrom}.locus --chrom 2 --chr_out $chrom \
--pos 3 --ident 1 --ref 4 --alt 5 --beta 8 --se 9 --indels N
qsub -t $chrom $scripts/sub_vcf_no_mono_1k.sh
}
get_chrom_ref 2 
get_chrom_ref 11 

function get_interval_ref 
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
Rscript --vanilla $scripts/extract_snp_interval.R chr${chrom}.locus $snp $interval
python $utils/create_batch_scripts.py $scripts/sub_CalcLD_1KG_VC_1K.sh  \
sub_CalcLD_${snp}_${interval}.sh my_file=$analysis/chr${chrom}.${snp}.${interval}.locus chrom=$chrom
qsub sub_CalcLD_${snp}_${interval}.sh
}

for chrom in {2..22}
do
python $scripts/tab2locus_paintor.py --tab $gwas/results.euro.pval.1k \
--out $analysis/chr${chrom}.locus2 --chrom 2 --chr_out $chrom \
--pos 3 --ident 12 --ref 4 --alt 5 --beta 8 --se 9 --indels Y
done


for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 5000 50000 250000 500000 1500000 
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
Rscript --vanilla $scripts/extract_snp_interval.R chr${chrom}.locus2 $snp $interval
python $utils/create_batch_scripts.py $scripts/sub_CalcLD_1KG_VC_1K.sh  \
sub_CalcLD_${snp}_${interval}.sh my_file=$analysis/chr${chrom}.${snp}.${interval}.locus2 chrom=$chrom
qsub sub_CalcLD_${snp}_${interval}.sh
done
done


#calculate LD for all SNPs and intervals considered.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do 
for my_int in 5000 50000 250000 500000 1500000 
do 
qsub -v reference=$gwas/results.euro.pval.1k,snp=$my_snps,interval=$my_int $utils/sub_extract_snp_interval.sh
done


#Carry out analysis using haploblocks generated by BigLD on 1K EUR data.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
qsub -v my_rsid=$rsid,my_chrom=$chrom,my_pos=$pos,my_bigld=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bigld,my_gwas=$gwas/results.euro.pval.1k $scripts_ref/sub_extract_snp_interval_bigld.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt


get_interval_ref rs61813875 1500000
get_interval_ref rs112111458 1500000
get_interval_ref rs2212434 1500000

#Going to use the following additional intervals from T. Bettram's analysis:
#500kbp
#280kbp
#100kbp - this one should capture the majority of variants
get_interval_ref rs61813875 250000
get_interval_ref rs112111458 250000
get_interval_ref rs2212434 250000

get_interval_ref rs61813875 140000
get_interval_ref rs112111458 140000
get_interval_ref s2212434 140000

get_interval_ref rs61813875 50000
get_interval_ref rs112111458 50000
get_interval_ref rs2212434 50000


#Get all the remaining chromosomes.
function get_chrom_ref 
{
chrom=$1
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
qsub -t $chrom $scripts/sub_vcf_no_mono_1k.sh
}

#Convert all chromosome to plink transposed format
for chrom in {1..22}
do
qsub -v input_file=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz,optional="--recode transpose" $utils/sub_vcf_to_plink.sh
done

#Convert all chromosomes to BIM file.
#Conversion removes padding 0s for chromosome numbers, for comparison with Eczema GWAS. 
for chrom in {1..22}
do
zcat $analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz \
| grep -v "#" | awk -v OFS="\t" '{print $1, $3, '0', $2, $4, $5}' | sed 's/^0//' >$analysis/diagnostics/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim
done

#Concatenate all bim files into 1.
touch $analysis/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim
for chrom in {1..22}
do
cat $analysis/diagnostics/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim >>$analysis/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim
done

#Marker QC diagnosis between GWAS and SNP datasets (both with 23&me and without it).
cd $analysis/diagnostics
for my_bim in *.bim
do
python $utils/compare_bim_files_pos.py $gwas/results.published.bim $my_bim >results.published_vs_${my_bim%.bim}.log
done
cd $analysis

#Generate GWAS table fixed in the following way:
#Remove variants where we have different alleles at the same positions, identified by the same rsids.
#Change the rsids of variants to match that in the 1k dataset in cases where we have the same position, same chromosome, same alleles but different rsids.
#Remove variants present in GWAS file but not in the reference.
#Following this procedure we can match the new table using rsids in any other datasets imputed to 1k.

python $utils/harmonize_rsid.py $analysis/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim $gwas/results.euro.bim >$gwas/results.euro.1k.bim
python $utils/update_rsid.py --bim $gwas/results.euro.1k.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$gwas/results.euro.pval.1k

#Harmonize the effect sizes so that they are relevant to the same reference allele as in the UKBiobank file. 
python $utils/harmonize_beta.py --tab $gwas/results.euro.pval.1k --ref $analysis/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim --header_tab Y --header_ref N \
--rsid_tab 12 --rsid_ref 2 --effect_tab 4 --alt_tab 5 --effect_ref 5 --beta_tab 8 --zscore_tab 10 --out $gwas/results.euro.pval.1k.harmonized_beta


#Generate GWAS table fixed in the following way:
#Keep all the variants and update the ID column with rsids of those found to be matching in the ref panel.
python $utils/update_rsid_all.py --bim $analysis/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 --rsid 1 >$gwas/results.euro.pval.all.1k

#Calculate allele frequencies
for chrom in {1..22}
do
vcftools --gzvcf $analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz --freq --out 1k.EUR.${chrom}
done

#Generate haploblocks using the European reference dataset, without indels.
for chrom in {1..22}
do
qsub -v file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz $scripts/sub_generate_BigLD_input.sh
done

for chrom in {1..22}
do
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono $scripts/sub_bigld.sh
done

qsub -v input_file=test,heuristic="nonhrst" $scripts/sub_bigld_gpart.sh

#For some reason, despite regular filesize, and despite increased memory up to 5-fold, runs on chromosomes 6 and 16 keep crashing.
#Going to create a separate run, with smaller file size - cut down in relation to our index SNP so that haploblock definition in the region should be correct.
rsid=rs4713555
chrom=6
pos=32575524
#Chromosome 6: range: 20000000:40000000
vcftools --gzvcf ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz --chr 6 --from-bp 20000000 --to-bp 40000000 --recode --stdout | gzip -c > ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp.vcf.gz

rsid=rs2041733
chrom=16
pos=11229589
#Chromosome 16: range: 1:20000000
vcftools --gzvcf ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz --chr 16 --from-bp 1 --to-bp 20000000 --recode --stdout | gzip -c > ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf_20mbp.gz 

qsub -v file=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp.vcf.gz $scripts/sub_generate_BigLD_input.sh
qsub -v file=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp.vcf.gz $scripts/sub_generate_BigLD_input.sh


qsub -v input_file=ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp $scripts/sub_bigld.sh
qsub -v input_file=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp $scripts/sub_bigld.sh

mv ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_20mbp.bigld ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bigld

#Run for chromosome got killed again - giving up.
#Contacted the author and they have a new gpart package, which allows a faster, heuristic-based estimation of LD blocks.
#Try it on our chromosomes - 6 and 16 which show problems with LD block generation.
rsid=rs4713555
chrom=6
pos=32575524
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,heuristic="fast" $scripts/sub_bigld_gpart.sh 

rsid=rs2041733
chrom=16
pos=11229589
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,heuristic="near-nonhrst" $scripts/sub_bigld_gpart.sh 

#Generated successful output for both chrom 6 and 16. Try chrom 6 with the near-nonheuristic algorithm.
rsid=rs4713555
chrom=6
pos=32575524
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,heuristic="near-nonhrst" $scripts/sub_bigld_gpart.sh 

#Haploblock generation using gpart and D' as metric.
for chrom in 1 2 3 4 5 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22
do
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,heuristic="nonhrst" $scripts/sub_bigld_gpart.sh 
done
#Two chromosomes (2 and 17) took too long to generate output. Run them again for longer.

#Generate interval from gpart, bigld partitioning results for finemap and paintor.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
qsub -v my_rsid=$rsid,my_chrom=$chrom,my_pos=$pos,my_bigld=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.gpart,my_gwas=chr${chrom}.locus2 $scripts/sub_extract_snp_interval_bigld_finemap_paintor.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt

while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
qsub -v my_rsid=$rsid,my_chrom=$chrom,my_pos=$pos,my_bigld=$analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bigld,my_gwas=chr${chrom}.locus2 $scripts/sub_extract_snp_interval_bigld_finemap_paintor2.sh
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Generate BED files with BigLD intervals.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
Rscript --vanilla $scripts/extract_bigld_bed.R $rsid $chrom $pos $analysis/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bigld
done < $gwas/paternoster_2015_index_snps_sorted.txt


#Extract SNP interval for analysis using haploblocks from gpart and bigld.
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in gpart bigld
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
python $utils/create_batch_scripts.py $scripts/sub_CalcLD_1KG_VC_1K.sh  \
sub_CalcLD_${snp}_${interval}.sh my_file=$analysis/${snp}.${interval}.locus2 chrom=$chrom
qsub sub_CalcLD_${snp}_${interval}.sh
done
done

#In certain cases, SNPs of interest have not be assigned to any haploblock. Plot their LD pattern.
for snp in rs112111458 rs12295535 rs10791824 rs4643526 rs61813875 rs7127307 rs77714197
do
for interval in 5000 50000 
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
pos=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f3)
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,my_results=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.gpart,my_chrom=$chrom,snp_id=$snp,my_position=$pos,my_interval=$interval,my_dataset="1k_EUR" $scripts/sub_bigld_gpart_plotld.sh 
done
done

#Better yet, plot LD patterns of all SNPs in the analysis
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 5000 50000 
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
pos=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f3)
qsub -v input_file=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono,my_results=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.gpart,my_chrom=$chrom,snp_id=$snp,my_position=$pos,my_interval=$interval,my_dataset="1k_EUR" $scripts/sub_bigld_gpart_plotld.sh 
done
done

#Create a subset of the VCF files within 1500 kbp around the focal SNP.
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
my_start=$(expr $pos - 1500000)
my_end=$(expr $pos + 1500000)
echo $chrom
echo $my_start
echo $my_end
vcftools --gzvcf ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz --chr $chrom --from-bp $my_start --to-bp $my_end --recode --stdout | gzip -c > ${rsid}_3Mbp_1kEUR.vcf.gz
done < $gwas/paternoster_2015_index_snps_sorted.txt

#Calculate LD based on the files.
for my_snps in *_3Mbp_1kEUR.vcf.gz
do
qsub -v my_file=$my_snps $utils/sub_vcf_ld.sh 
done

#Convert the Ref Panel to bed/bim/fam Plink format to be used in regfm pipeline.
for a in ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz
do
qsub -v input_file=$a $utils/sub_vcf_to_plink.sh
done

#Move the original VCF files and filtered for monomorphic in EUR to the new folder for the MendelVar project
mv ALL.chr*EUR_no_mono.vcf.gz /newhome/qh18484/MendelVar/ref_panel
mv ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz /newhome/qh18484/MendelVar/ref_panel
mv integrated_call_samples_v3.20130502.ALL.panel /newhome/qh18484/MendelVar/ref_panel
