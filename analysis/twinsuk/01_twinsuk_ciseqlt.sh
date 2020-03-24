#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/twinsuk
scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk
skin_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/skin
lcl_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/LCLs
var_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/gen
sample_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/sample
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils

#Identification of cis-eQTLs (up to 1Mbp upstream of TSS) from TwinsUK for our 44 GWAS loci, for skin and LCLs.
#Consider all SNPs within 1Mbp of our index SNP.
cd $analysis

#Generate a TwinsUK version of the GWAS results. First: a BIM file.
for a in $var_data/*.gen
do
chrom=$(echo $a | cut -d"." -f2 | cut -d"r" -f2 | sed 's/^0//')
cat $a | awk -v OFS="\t" -v chr="$chrom" '{print chr, $2, '0', $3, $4, $5}' >>twinsuk.bim
done

python $utils/update_rsid.py --bim twinsuk.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$gwas/results.euro.pval.twinsuk

#Harmonize the effect sizes so that they are relevant to the same reference allele as in the TwinsUK file. 
python $utils/harmonize_beta.py --tab $gwas/results.euro.pval.twinsuk --ref twinsuk.bim  --header_tab Y --header_ref N \
--rsid_tab 12 --rsid_ref 2 --effect_tab 4 --alt_tab 5 --effect_ref 5 --beta_tab 8 --zscore_tab 10 --out $gwas/results.euro.pval.twinsuk.harmonized_beta

#Create BED files with intervals to select for each index SNP.
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-50000, $3+50000, $1}' >paternoster_2015_index_snps_sorted_100kbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-125000, $3+125000, $1}' >paternoster_2015_index_snps_sorted_250kbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-500000, $3+500000, $1}' >paternoster_2015_index_snps_sorted_1Mbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-250000, $3+250000, $1}' >paternoster_2015_index_snps_sorted_500kbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-1500000, $3+1500000, $1}' >paternoster_2015_index_snps_sorted_3Mbp.bed


#Convert bim to bed.
cat twinsuk.bim | awk -v OFS="\t" '{print $1, $4, $4, $2}' >twinsuk.bed
S
while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_3Mbp.bed
done < paternoster_2015_index_snps_sorted_3Mbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_1Mbp.bed
done < paternoster_2015_index_snps_sorted_1Mbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_500kbp.bed
done < paternoster_2015_index_snps_sorted_500kbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_100kbp.bed
done < paternoster_2015_index_snps_sorted_100kbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_250kbp.bed
done < paternoster_2015_index_snps_sorted_250kbp.bed

#Identify SNPs overlapping the interval
for a in rs*.bed
do
intersectBed -wa -a twinsuk.bed -b $a>${a%.bed}_snps.bed
done

#Identify genes overlapping the interval around each index SNP.
for a in rs*snps.bed
do
intersectBed -wa -a $HOME/analysis/colocalization/coloc/sun_pqtl/gencode.v19.annotation.gtf -b $a >${a%.bed}.gtf
done

#Get a list of ENG Ids to be extraced in each interval.
for a in rs*_snps.gtf
do
cat $a | sed 's/^.*\(gene_id .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | cut -d"." -f1 | sort | uniq >${a%.gtf}.gene_ids
done

#Going to carry out analysis on the 100kbp, 250kbp and 1Mbp and 3Mbp interval for now. 
#extract transcript data.
for b in rs*1Mbp_snps.gene_ids rs*3Mbp_snps.gene_ids
do
while read line
do
set -- $linels
gene_id=$1
grep -w "$gene_id" $skin_data/data.rpkm >>skin/${b%.gene_ids}.rpkm
done < $b
done

for b in rs*1Mbp_snps.gene_ids rs*3Mbp_snps.gene_ids
do
while read line
do
set -- $line
gene_id=$1
grep -w "$gene_id" $lcl_data/data.rpkm >>lcl/${b%.gene_ids}.rpkm
done < $b

#Fix header in each file.
for a in skin/*_snps.rpkm
do
head -1 $skin_data/data.rpkm | cat - $a >${a%.rpkm}_headers.rpkm
done

for a in lcl/*_snps.rpkm
do
head -1 $lcl_data/data.rpkm | cat - $a >${a%.rpkm}_headers.rpkm
done

#Subset gene expression results to asthma and non-asthma only individuals.
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cd $analysis
cd skin
Rscript --vanilla $scripts/subset_eczema.R ${snp}_${interval}_snps_headers.rpkm $analysis/chr2_Eurobats_Public_alt_pheno.sample
cd $analysis
cd lcl
Rscript --vanilla $scripts/subset_eczema.R ${snp}_${interval}_snps_headers.rpkm $analysis/chr2_Eurobats_Public_alt_pheno.sample
done
done

#Create files with list of individuals in each subset: eczema, noneczema, all.
for a in skin lcl
do
cd $analysis/$a 
head -1 rs112111458_3Mbp_snps_headers_all.rpkm | tr '\t' '\n' | tail -n +2 >expression_individuals_all
head -1 rs112111458_3Mbp_snps_headers_all_eczema.rpkm | tr '\t' '\n' | tail -n +2 >expression_individuals_eczema
head -1 rs112111458_3Mbp_snps_headers_all_noneczema.rpkm | tr '\t' '\n' | tail -n +2 >expression_individuals_noneczema
done

#Generate same files, but to be used by plink.
for a in skin lcl
do
cd $analysis/$a 
awk -v OFS="\t" '{print $1, $1 }' expression_individuals_all >expression_individuals_all_plink
awk -v OFS="\t" '{print $1, $1 }' expression_individuals_eczema >expression_individuals_eczema_plink
awk -v OFS="\t" '{print $1, $1 }' expression_individuals_noneczema >expression_individuals_noneczema_plink
done

#Create files with SNP ids in each interval.
for my_file in rs*_snps.bed
do
cut -f4 $my_file > ${my_file%.bed}.lst
done

#Subset genotype file.
cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
if [ "$chrom" -ge 1 -a "$chrom" -le 9 ]
then
	my_chrom=$(echo 0$chrom)
else 
	my_chrom=$chrom
fi
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/skin/${snp}_${interval}_all.gen \
--os $analysis/skin/${snp}_${interval}_all.sample --sample_id $analysis/skin/expression_individuals_all --inclusion $analysis/${snp}_${interval}_snps.lst
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/skin/${snp}_${interval}_eczema.gen \
--os $analysis/skin/${snp}_${interval}_eczema.sample --sample_id $analysis/skin/expression_individuals_eczema --inclusion $analysis/${snp}_${interval}_snps.lst
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/skin/${snp}_${interval}_noneczema.gen \
--os $analysis/skin/${snp}_${interval}_noneczema.sample --sample_id $analysis/skin/expression_individuals_noneczema --inclusion $analysis/${snp}_${interval}_snps.lst
done
done

cd $analysis/lcl 
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
if [ "$chrom" -ge 1 -a "$chrom" -le 9 ]
then
	my_chrom=$(echo 0$chrom)
else 
	my_chrom=$chrom
fi
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/lcl/${snp}_${interval}_all.gen \
--os $analysis/lcl/${snp}_${interval}_all.sample --sample_id $analysis/lcl/expression_individuals_all --inclusion $analysis/${snp}_${interval}_snps.lst
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/lcl/${snp}_${interval}_eczema.gen \
--os $analysis/lcl/${snp}_${interval}_eczema.sample --sample_id $analysis/lcl/expression_individuals_eczema --inclusion $analysis/${snp}_${interval}_snps.lst
gtool -S --g $var_data/data.chr$my_chrom.gen --s $sample_data/data.chr$my_chrom.sample --og $analysis/lcl/${snp}_${interval}_noneczema.gen \
--os $analysis/lcl/${snp}_${interval}_noneczema.sample --sample_id $analysis/lcl/expression_individuals_noneczema --inclusion $analysis/${snp}_${interval}_snps.lst
done
done

#Merge all gen files into 1.
for a in $var_data/data.*.gen
do
cat $a >> data.all.gen
done

#Generate PCA matrix of genetic relationship with SNPRelate.
#Convert chromosome 1 to Plink format.
plink --gen $var_data/data.chr01.gen --sample $sample_data/data.chr01.sample --oxford-single-chr 1 --make-bed --out data.chr01
#Convert sll chromosomes to Plink format.
plink --gen data.all.gen --sample $sample_data/data.chr01.sample --oxford-single-chr 1 --make-bed --out data.all
qsub $scripts/sub_twinsuk_to_plink.sh
#Calculate PCA
qsub $scripts/sub_pca_twinsuk.sh
qsub $scripts/sub_pca_twinsuk_chr1.sh

#Create covariate file with sex, age and genetic covariates (6 top PCs) for all individuals, eczema and noneczema 
#for matrixeqtl and PEER.
cd $analysis/skin
Rscript --vanilla $scripts/matrix_eqtl_covariate.R
cd $analysis/lcl
Rscript --vanilla $scripts/matrix_eqtl_covariate.R

#Get files with individuals to be analyzed in the matrixeqtl analysis (both gene expression and genotypes)
cd $analysis/skin
for a in all eczema noneczema
do
head -1 twinsuk_matrixeqtl_covariates.${a} | tr '\t' '\n' >twinsuk_individuals_${a}
done

cd $analysis/lcl
for a in all eczema noneczema
do
head -1 twinsuk_matrixeqtl_covariates.${a} | tr '\t' '\n' >twinsuk_individuals_${a}
done

#Generate kinship matrix.
cd $analysis/skin
for a in all eczema noneczema
do
Rscript --vanilla $scripts/generate_kinship_matrix.R twinsuk_individuals_${a} twinsuk_relatedness_${a} 
done

cd $analysis/lcl
for a in all eczema noneczema
do
Rscript --vanilla $scripts/generate_kinship_matrix.R twinsuk_individuals_${a} twinsuk_relatedness_${a} 
done

#Convert genotypes and expression values to MatrixeQTL format.
#Note added 04/01/2018 - the following conversion to MatrixeQTL format is incorrect, hence the final results below need to be ignored (as I did earlier anyway). The problem lies in not recognizing that the gen/sample format provides probabilities of each of the three genotypes per individual, while the MatrixeQTL format looks at the dosages of the reference allele per individual, i.e. one value per individual.
#Normalize RPKM values to standard normal according to matrixeqtl, for use in Peer. 
#First, add missing tab in the beginning of the header.
cd $analysis/skin
#sed -i '1s/^/\t/' $skin_data/data.rpkm 
qsub -v input=$skin_data/data.rpkm $scripts/sub_normalize_expression.sh 

cd $analysis/lcl
#sed -i '1s/^/\t/' $lcl_data/data.rpkm 
qsub -v input=$lcl_data/data.rpkm $scripts/sub_normalize_expression.sh 


cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_all.gen ${snp}_${interval}_all.sample \
${snp}_${interval}_snps_headers_all.rpkm twinsuk_individuals_all
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_eczema.gen ${snp}_${interval}_eczema.sample \
${snp}_${interval}_snps_headers_eczema.rpkm twinsuk_individuals_eczema
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_noneczema.gen ${snp}_${interval}_noneczema.sample \
${snp}_${interval}_snps_headers_noneczema.rpkm twinsuk_individuals_noneczema
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_all.gen ${snp}_${interval}_all.sample \
${snp}_${interval}_snps_headers_all.rpkm twinsuk_individuals_all
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_eczema.gen ${snp}_${interval}_eczema.sample \
${snp}_${interval}_snps_headers_eczema.rpkm twinsuk_individuals_eczema
Rscript --vanilla $scripts/matrix_eqtl_input.R ${snp}_${interval}_noneczema.gen ${snp}_${interval}_noneczema.sample \
${snp}_${interval}_snps_headers_noneczema.rpkm twinsuk_individuals_noneczema
done
done

#Remove Sex from covariates as it is collinear - only females in the dataset
#Run MatrixEQTL: using only age as covariate and relatednesss as error matrix.
cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression twinsuk_matrixeqtl_covariates_nopca.all \
${snp}_${interval}_all_rel.matrixeqtl twinsuk_relatedness_all
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
twinsuk_matrixeqtl_covariates_nopca.eczema ${snp}_${interval}_eczema_rel.matrixeqtl twinsuk_relatedness_eczema
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
twinsuk_matrixeqtl_covariates_nopca.noneczema ${snp}_${interval}_noneczema_rel.matrixeqtl twinsuk_relatedness_noneczema
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression twinsuk_matrixeqtl_covariates_nopca.all \
${snp}_${interval}_all_rel.matrixeqtl twinsuk_relatedness_all
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
twinsuk_matrixeqtl_covariates_nopca.eczema ${snp}_${interval}_eczema_rel.matrixeqtl twinsuk_relatedness_eczema
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
twinsuk_matrixeqtl_covariates_nopca.noneczema ${snp}_${interval}_noneczema_rel.matrixeqtl twinsuk_relatedness_noneczema
done
done

#Run Peer analysis
cd $analysis/skin
qsub -v expression=data_all_normal.rpkm,covariate=twinsuk_peer_covariates.all $scripts/sub_run_peer.sh 
qsub -v expression=data_eczema_normal.rpkm,covariate=twinsuk_peer_covariates.eczema $scripts/sub_run_peer.sh
qsub -v expression=data_noneczema_normal.rpkm,covariate=twinsuk_peer_covariates.noneczema $scripts/sub_run_peer.sh


cd $analysis/lcl
qsub -v expression=data_all_normal.rpkm,covariate=twinsuk_peer_covariates.all $scripts/sub_run_peer.sh 
qsub -v expression=data_eczema_normal.rpkm,covariate=twinsuk_peer_covariates.eczema $scripts/sub_run_peer.sh
qsub -v expression=data_noneczema_normal.rpkm,covariate=twinsuk_peer_covariates.noneczema $scripts/sub_run_peer.sh


#Create covariate file with results including hidden confounders from Peer and known covariates.
cd $analysis/skin 
Rscript --vanilla $scripts/prepare_peer_covariates_disease.R twinsuk_matrixeqtl_covariates.noneczema test_noneczema_normal_hidden_factors
Rscript --vanilla $scripts/prepare_peer_covariates_disease.R twinsuk_matrixeqtl_covariates.eczema test_eczema_normal_hidden_factors
Rscript --vanilla $scripts/prepare_peer_covariates_all.R twinsuk_matrixeqtl_covariates.all test_all_normal_hidden_factors

cd $analysis/lcl 
Rscript --vanilla $scripts/prepare_peer_covariates_disease.R twinsuk_matrixeqtl_covariates.noneczema test_noneczema_normal_hidden_factors
Rscript --vanilla $scripts/prepare_peer_covariates_disease.R twinsuk_matrixeqtl_covariates.eczema test_eczema_normal_hidden_factors
Rscript --vanilla $scripts/prepare_peer_covariates_all.R twinsuk_matrixeqtl_covariates.all test_all_normal_hidden_factors

#Run MatrixEQTL: using age, disease status, PCs along with Peer factors as covariates.
cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
#Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression test_all_normal_hidden_factors_pc_covariates \
#${snp}_${interval}_all_pc_peer.matrixeqtl 
Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
test_eczema_normal_hidden_factors_pc_covariates ${snp}_${interval}_eczema_pc_peer.matrixeqtl 
Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
test_noneczema_normal_hidden_factors_pc_covariates ${snp}_${interval}_noneczema_pc_peer.matrixeqtl 
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
#Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression test_all_normal_hidden_factors_pc_covariates \
#${snp}_${interval}_all_pc_peer.matrixeqtl 
Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
test_eczema_normal_hidden_factors_pc_covariates ${snp}_${interval}_eczema_pc_peer.matrixeqtl 
Rscript --vanilla $scripts/matrix_eqtl.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
test_noneczema_normal_hidden_factors_pc_covariates ${snp}_${interval}_noneczema_pc_peer.matrixeqtl 
done
done


#Run MatrixEQTL: Using age, disease, Peer covariates as well as relatedness matrix.
cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
#Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression test_all_normal_hidden_factors_nopc_covariates \
#${snp}_${interval}_all_pc_rel.matrixeqtl twinsuk_relatedness_all
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
test_eczema_normal_hidden_factors_nopc_covariates ${snp}_${interval}_eczema_pc_rel.matrixeqtl twinsuk_relatedness_eczema
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
test_noneczema_normal_hidden_factors_nopc_covariates ${snp}_${interval}_noneczema_pc_rel.matrixeqtl twinsuk_relatedness_noneczema
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
#Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_all.genotype ${snp}_${interval}_all.expression test_all_normal_hidden_factors_nopc_covariates \
#${snp}_${interval}_all_pc_rel.matrixeqtl twinsuk_relatedness_all
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_eczema.genotype ${snp}_${interval}_eczema.expression \
test_eczema_normal_hidden_factors_nopc_covariates ${snp}_${interval}_eczema_pc_rel.matrixeqtl twinsuk_relatedness_eczema
Rscript --vanilla $scripts/matrix_eqtl_relate.R ${snp}_${interval}_noneczema.genotype ${snp}_${interval}_noneczema.expression \
test_noneczema_normal_hidden_factors_nopc_covariates ${snp}_${interval}_noneczema_pc_rel.matrixeqtl twinsuk_relatedness_noneczema
done
done


#Check which protocol: pc_peer, rel_peer or rel results in optimal number of significant eQTLs.
cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cat ${snp}_${interval}_eczema_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_rel.matrixeqtl.sig
cat ${snp}_${interval}_noneczema_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>noneczema_rel.matrixeqtl.sig 
cat ${snp}_${interval}_all_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_rel.matrixeqtl.sig
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
#cat ${snp}_${interval}_eczema_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_rel.matrixeqtl.sig
#cat ${snp}_${interval}_noneczema_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.15) print $0}' >>noneczema_rel.matrixeqtl.sig 
cat ${snp}_${interval}_all_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_rel.matrixeqtl.sig
done
done

cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cat ${snp}_${interval}_eczema_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_pc_peer.matrixeqtl.sig
cat ${snp}_${interval}_noneczema_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>noneczema_pc_peer.matrixeqtl.sig 
#cat ${snp}_${interval}_all_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_pc_peer.matrixeqtl.sig
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cat ${snp}_${interval}_eczema_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_pc_peer.matrixeqtl.sig
cat ${snp}_${interval}_noneczema_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>noneczema_pc_peer.matrixeqtl.sig 
#cat ${snp}_${interval}_all_pc_peer.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_pc_peer.matrixeqtl.sig
done
done


cd $analysis/skin
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cat ${snp}_${interval}_eczema_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_pc_rel.matrixeqtl.sig
cat ${snp}_${interval}_noneczema_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>noneczema_pc_rel.matrixeqtl.sig 
#cat ${snp}_${interval}_all_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_pc_rel.matrixeqtl.sig
done
done

cd $analysis/lcl
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
cat ${snp}_${interval}_eczema_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>eczema_pc_rel.matrixeqtl.sig
cat ${snp}_${interval}_noneczema_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>noneczema_pc_rel.matrixeqtl.sig 
#cat ${snp}_${interval}_all_pc_rel.matrixeqtl.matrixeqtl | awk -F '\t' '{if ($5 < 0.05) print $0}' >>all_pc_rel.matrixeqtl.sig
done
done

#Tally the numbers
#For skin
for a in rel pc_peer pc_rel #A more suitable name for pc_rel would be peer_rel, as we are using Peer factors and relatedness matrix here.
do
echo $a
echo "eczema"
cat eczema_${a}.matrixeqtl.sig | cut -f1,2 | sort | uniq | wc -l
echo "noneczema"
cat noneczema_${a}.matrixeqtl.sig | cut -f1,2 | sort | uniq | wc -l
done

rel
eczema
1127
noneczema
1470

pc_peer
eczema
322
noneczema
345

pc_rel
eczema
257
noneczema
957

#In short, the rel protocol with no Peer correction is superior. 
#Check if we get more looking at a dataset as a whole.
cat all_rel.matrixeqtl.sig | cut -f1,2 | sort | uniq | wc -l
#Surprisingly fewer: just 601!

#For LCLs.
rel
eczema
102
noneczema
999
pc_peer
eczema
0
noneczema
263
pc_rel
eczema
148
noneczema
1121

#Not such a big difference here, rel and pc_rel protocols comparable
#so sticking with rel.

cat all_rel.matrixeqtl.sig | cut -f1,2 | sort | uniq | wc -l
#Surprisingly, here we have more when looking at combined all file.
#1270

