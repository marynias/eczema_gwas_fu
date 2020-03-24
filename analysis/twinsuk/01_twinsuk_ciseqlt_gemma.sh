#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/twinsuk
scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk
coloc_scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
skin_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/skin
lcl_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/LCLs
var_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/gen
sample_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/sample
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils

#Identification of cis-eQTLs (up to 1Mbp upstream of TSS) from TwinsUK for our 44 GWAS loci, for skin and LCLs.
#Consider all SNPs within 1Mbp of our index SNP. Using gemma to idenitfy eQTLs here.
#Just skin analysis here.
cd $analysis/gemma/skin
#Convert gen files to bimbam
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864 rs112111458
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
my_samples=$(wc -l $analysis/skin/${snp}_${interval}_all.sample | cut -d" " -f1)
sample_no="$((my_samples-2))"
qsub -v gen_file=$analysis/skin/${snp}_${interval}_all.gen,number_of_samples=$sample_no,bb_file=${snp}_${interval}_all.bimbam $utils/sub_gen2bimbam.sh
my_samples=$(wc -l $analysis/skin/${snp}_${interval}_eczema.sample | cut -d" " -f1)
sample_no="$((my_samples-2))"
qsub -v gen_file=$analysis/skin/${snp}_${interval}_eczema.gen,number_of_samples=$sample_no,bb_file=${snp}_${interval}_eczema.bimbam $utils/sub_gen2bimbam.sh
my_samples=$(wc -l $analysis/skin/${snp}_${interval}_noneczema.sample | cut -d" " -f1)
sample_no="$((my_samples-2))"
qsub -v gen_file=$analysis/skin/${snp}_${interval}_noneczema.gen,number_of_samples=$sample_no,bb_file=${snp}_${interval}_noneczema.bimbam $utils/sub_gen2bimbam.sh
done
done

#Prepare a list of order of samples in each subset
tail -n +3 $analysis/skin/rs2227483_3Mbp_all.sample | cut -d" " -f1 >genosample_all
tail -n +3 $analysis/skin/rs2227483_3Mbp_eczema.sample | cut -d" " -f1 >genosample_eczema
tail -n +3 $analysis/skin/rs2227483_3Mbp_noneczema.sample | cut -d" " -f1 >genosample_noneczema


#Transform the values with rnorm function in GenABEL.

#Create relatedness matrix.
#Subset all chromosome data to all, eczema and non-eczema samples for skin
for i in {01..22}
do
gtool -S --g $var_data/data.chr${i}.gen \
--s $sample_data/data.chr${i}.sample \
--og $analysis/gemma/skin/skin_all.chr${i}.gen \
--os $analysis/gemma/skin/skin_all.chr${i}.sample \
--sample_id $analysis/gemma/skin/genosample_all \
--log  $analysis/gemma/skin/skin_all.chr${i}.log 
done


for i in {01..22}
do
gtool -S --g $var_data/data.chr${i}.gen \
--s $sample_data/data.chr${i}.sample \
--og $analysis/gemma/skin/skin_eczema.chr${i}.gen \
--os $analysis/gemma/skin/skin_eczema.chr${i}.sample \
--sample_id $analysis/gemma/skin/genosample_eczema \
--log  $analysis/gemma/skin/skin_eczema.chr${i}.log 
done

for i in {01..22}
do
gtool -S --g $var_data/data.chr${i}.gen \
--s $sample_data/data.chr${i}.sample \
--og $analysis/gemma/skin/skin_noneczema.chr${i}.gen \
--os $analysis/gemma/skin/skin_noneczema.chr${i}.sample \
--sample_id $analysis/gemma/skin/genosample_noneczema \
--log  $analysis/gemma/skin/skin_noneczema.chr${i}.log 
done


#convert imputation .gen file to BIMBAM format (mean_genotype file).
for i in {01..22}
do
qsub -v number_of_samples=672,gen_file=$analysis/gemma/skin/skin_all.chr${i}.gen,bb_file=$analysis/gemma/skin/skin_all.chr${i}.bimbam $utils/sub_gen2bimbam.sh
done

for i in {01..22}
do
qsub -v number_of_samples=86,gen_file=$analysis/gemma/skin/skin_eczema.chr${i}.gen,bb_file=$analysis/gemma/skin/skin_eczema.chr${i}.bimbam $utils/sub_gen2bimbam.sh
done

for i in {01..22}
do
qsub -v number_of_samples=530,gen_file=$analysis/gemma/skin/skin_noneczema.chr${i}.gen,bb_file=$analysis/gemma/skin/skin_noneczema.chr${i}.bimbam $utils/sub_gen2bimbam.sh
done

#merge into one file
cat $analysis/gemma/skin/skin_all.chr*.bimbam > all.bimbam 
cat $analysis/gemma/skin/skin_noneczema.chr*.bimbam > noneczema.bimbam 
cat $analysis/gemma/skin/skin_eczema.chr*.bimbam > eczema.bimbam 

#create fake pheno_file where there is no missing phenotype data (else those samples will not be included when computing the relatedness matrix)
for i in $(seq 1 672) ; do echo 1; done > all.pheno
for i in $(seq 1 86) ; do echo 1; done > eczema.pheno
for i in $(seq 1 530) ; do echo 1; done > noneczema.pheno

#create relatedness matrix. Warning! Gemma puts all the results into newly created output directory,
qsub $scripts/gemma_all_rel.sh
qsub $scripts/gemma_eczema_rel.sh
qsub $scripts/gemma_noneczema_rel.sh

#Prepare SNP annotation file (SNP id, SNP position and SNP chrom)
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
python $scripts/gemma_snp_annotation.py $analysis/twinsuk.bim $analysis/gemma/skin/${snp}_${interval}_all.bimbam
python $scripts/gemma_snp_annotation.py $analysis/twinsuk.bim $analysis/gemma/skin/${snp}_${interval}_eczema.bimbam
python $scripts/gemma_snp_annotation.py $analysis/twinsuk.bim $analysis/gemma/skin/${snp}_${interval}_noneczema.bimbam
done
done

#Prepare gene expression file, using files generated in the $analysis/skin folder.
Rscript $scripts/genabel.R

#run eQTL analysis 
# "-n 2" will tell it to look at the second column of the pheno file that has the expression data
# use "-n 1 2" to look at both expression and case-control status
for my_file in *_eczema.rpkm.rnt  
do
qsub -v file=$my_file,factors=1 $scripts/sub_gemma_eqtl.sh
done

for my_file in *_noneczema.rpkm.rnt  
do
qsub -v file=$my_file,factors=1 $scripts/sub_gemma_eqtl.sh
done

#Additional factor - disease status for all samples
for my_file in *_*all.rpkm.rnt  
do
qsub -v file=$my_file,factors=2 $scripts/sub_gemma_eqtl.sh
done

for file in *_eczema.rpkm.rnt  
do
gene=$(echo $file | cut -d"_" -f1)
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
group=$(echo $file | cut -d"_" -f6 | cut -d"." -f1)
gemma -g ${snp}_${interval}_${group}.bimbam -p $file -n 1 -a ${snp}_${interval}_${group}.snps -c gemma_${group}_covariate.txt -lmm 4 -o ${gene}_${snp}_${interval}_${group}.gemma -k output/${group}_relatedness.cXX.txt
done

for file in *_noneczema.rpkm.rnt  
do
gene=$(echo $file | cut -d"_" -f1)
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
group=$(echo $file | cut -d"_" -f6 | cut -d"." -f1)
gemma -g ${snp}_${interval}_${group}.bimbam -p $file -n 1 -a ${snp}_${interval}_${group}.snps -c gemma_${group}_covariate.txt -lmm 4 -o ${gene}_${snp}_${interval}_${group}.gemma -k output/${group}_relatedness.cXX.txt
done

for file in *_all.rpkm.rnt 
do
gene=$(echo $file | cut -d"_" -f1)
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
group=$(echo $file | cut -d"_" -f6 | cut -d"." -f1)
gemma -g ${snp}_${interval}_${group}.bimbam -p $file -n 2 -a ${snp}_${interval}_${group}.snps -c gemma_${group}_covariate.txt -lmm 4 -o ${gene}_${snp}_${interval}_${group}.gemma -k output/${group}_relatedness.cXX.txt
done


#Run coloc.
for file in output/*all.gemma.assoc.txt
do
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
Rscript $scripts/twinsuk_gemma_coloc_input.R $analysis/rel/${snp}_${interval}.gwas $file 672
done

for file in output/*_eczema.gemma.assoc.txt
do
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
Rscript $scripts/twinsuk_gemma_coloc_input.R $analysis/rel/${snp}_${interval}.gwas $file 86
done

for file in output/*_noneczema.gemma.assoc.txt
do
snp=$(echo $file | cut -d"_" -f2)
interval=$(echo $file | cut -d"_" -f3)
Rscript $scripts/twinsuk_gemma_coloc_input.R $analysis/rel/${snp}_${interval}.gwas $file 530
done

#Repeat the same analysis for lcl, using the following number of samples in each category:
coloc_all lcl 764
coloc_eczema lcl 92
coloc_noneczema lcl 610

#Check the colocalisation results.
function pvals {
my_tissue=$1
cd $analysis/gemma/$my_tissue/output
for a in *totalb
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a
done
}

pvals skin
pvals lcl

#Skin
AP003774.1_rs10791824_3Mbp_noneczema.gemma.assoc.totalb
CD207_rs112111458_100kbp_all.gemma.assoc.totalb
CD207_rs112111458_100kbp_noneczema.gemma.assoc.totalb
CD207_rs112111458_1Mbp_all.gemma.assoc.totalb
CD207_rs112111458_1Mbp_noneczema.gemma.assoc.totalb
CD207_rs112111458_250kbp_all.gemma.assoc.totalb
CD207_rs112111458_250kbp_noneczema.gemma.assoc.totalb
CD207_rs112111458_3Mbp_all.gemma.assoc.totalb
CD207_rs112111458_3Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_100kbp_all.gemma.assoc.totalb
KIAA0391_rs2038255_100kbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_1Mbp_all.gemma.assoc.totalb
KIAA0391_rs2038255_1Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_250kbp_all.gemma.assoc.totalb
KIAA0391_rs2038255_250kbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_3Mbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_100kbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_100kbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2143950_1Mbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_1Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2143950_250kbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_250kbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2143950_3Mbp_all.gemma.assoc.totalb
KIAA1841_rs4643526_250kbp_all.gemma.assoc.totalb
LIME1_rs4809219_1Mbp_all.gemma.assoc.totalb
LIME1_rs4809219_1Mbp_noneczema.gemma.assoc.totalb
LIME1_rs4809219_250kbp_all.gemma.assoc.totalb
LIME1_rs4809219_250kbp_noneczema.gemma.assoc.totalb
LIME1_rs4809219_3Mbp_all.gemma.assoc.totalb
MRGBP_rs4809219_3Mbp_all.gemma.assoc.totalb
PRR5L_rs12295535_100kbp_all.gemma.assoc.totalb
PRR5L_rs12295535_100kbp_noneczema.gemma.assoc.totalb
PRR5L_rs12295535_1Mbp_all.gemma.assoc.totalb
PRR5L_rs12295535_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs12295535_250kbp_all.gemma.assoc.totalb
PRR5L_rs12295535_250kbp_noneczema.gemma.assoc.totalb
PRR5L_rs12295535_3Mbp_all.gemma.assoc.totalb
PRR5L_rs12295535_3Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_100kbp_all.gemma.assoc.totalb
PRR5L_rs2218565_100kbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_1Mbp_all.gemma.assoc.totalb
PRR5L_rs2218565_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_250kbp_all.gemma.assoc.totalb
PRR5L_rs2218565_250kbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_3Mbp_all.gemma.assoc.totalb
PRR5L_rs2218565_3Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_100kbp_all.gemma.assoc.totalb
PRR5L_rs2592555_100kbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_1Mbp_all.gemma.assoc.totalb
PRR5L_rs2592555_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_250kbp_all.gemma.assoc.totalb
PRR5L_rs2592555_250kbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_3Mbp_all.gemma.assoc.totalb
PRR5L_rs2592555_3Mbp_noneczema.gemma.assoc.totalb
RIIAD1_rs7512552_3Mbp_all.gemma.assoc.totalb
S100A12_rs61813875_3Mbp_all.gemma.assoc.totalb
S100A12_rs61813875_3Mbp_noneczema.gemma.assoc.totalb
SNX6_rs2038255_1Mbp_noneczema.gemma.assoc.totalb
SNX6_rs2038255_3Mbp_noneczema.gemma.assoc.totalb
SNX6_rs2143950_1Mbp_noneczema.gemma.assoc.totalb
SNX6_rs2143950_3Mbp_noneczema.gemma.assoc.totalb
STMN3_rs4809219_100kbp_all.gemma.assoc.totalb
STMN3_rs4809219_1Mbp_all.gemma.assoc.totalb
STMN3_rs4809219_250kbp_all.gemma.assoc.totalb
STMN3_rs4809219_3Mbp_all.gemma.assoc.totalb
TARS2_rs7512552_1Mbp_all.gemma.assoc.totalb
TARS2_rs7512552_3Mbp_all.gemma.assoc.totalb

while read line
do
set -- $line
my_file=$1
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
interval=$(echo $my_file | cut -d"_" -f3)
treatment=$(echo $my_file | cut -d"_" -f4 | cut -d"." -f1)
target_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${my_file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $my_file)
pval_p=$(awk 'NR == 7 {print $2}' ${my_file%b}p)
echo -e $snp'\t'$interval'\t'$treatment'\t'$gene'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>twinsuk_skin_coloc_results
done < significant_totalb.txt

#LCL
AQP10_rs12730935_3Mbp_all.gemma.assoc.totalb
BAG6_rs145809981_1Mbp_all.gemma.assoc.totalb
BAG6_rs145809981_1Mbp_noneczema.gemma.assoc.totalb
BAG6_rs41293864_1Mbp_all.gemma.assoc.totalb
BAG6_rs41293864_1Mbp_noneczema.gemma.assoc.totalb
CAVIN1_rs12951971_1Mbp_all.gemma.assoc.totalb
CAVIN1_rs12951971_1Mbp_noneczema.gemma.assoc.totalb
CAVIN1_rs17881320_1Mbp_all.gemma.assoc.totalb
CAVIN1_rs17881320_1Mbp_noneczema.gemma.assoc.totalb
ECM1_rs7512552_1Mbp_all.gemma.assoc.totalb
ECM1_rs7512552_1Mbp_noneczema.gemma.assoc.totalb
ECM1_rs7512552_3Mbp_all.gemma.assoc.totalb
ECM1_rs7512552_3Mbp_noneczema.gemma.assoc.totalb
FAM208B_rs6602364_1Mbp_all.gemma.assoc.totalb
FAM208B_rs6602364_1Mbp_noneczema.gemma.assoc.totalb
FAM208B_rs6602364_3Mbp_all.gemma.assoc.totalb
FAM208B_rs6602364_3Mbp_noneczema.gemma.assoc.totalb
GMEB2_rs4809219_1Mbp_all.gemma.assoc.totalb
GMEB2_rs4809219_1Mbp_noneczema.gemma.assoc.totalb
GMEB2_rs4809219_3Mbp_all.gemma.assoc.totalb
GMEB2_rs4809219_3Mbp_noneczema.gemma.assoc.totalb
IL18RAP_rs3917265_1Mbp_all.gemma.assoc.totalb
IL18RAP_rs3917265_1Mbp_noneczema.gemma.assoc.totalb
IL18RAP_rs3917265_3Mbp_all.gemma.assoc.totalb
IL18RAP_rs3917265_3Mbp_noneczema.gemma.assoc.totalb
IL18RAP_rs6419573_1Mbp_all.gemma.assoc.totalb
IL18RAP_rs6419573_1Mbp_noneczema.gemma.assoc.totalb
IL18RAP_rs6419573_3Mbp_all.gemma.assoc.totalb
IL18RAP_rs6419573_3Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_1Mbp_all.gemma.assoc.totalb
KIAA0391_rs2038255_1Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2038255_3Mbp_all.gemma.assoc.totalb
KIAA0391_rs2038255_3Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2143950_1Mbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_1Mbp_noneczema.gemma.assoc.totalb
KIAA0391_rs2143950_3Mbp_all.gemma.assoc.totalb
KIAA0391_rs2143950_3Mbp_noneczema.gemma.assoc.totalb
KRT20_rs12951971_3Mbp_all.gemma.assoc.totalb
KRT20_rs12951971_3Mbp_noneczema.gemma.assoc.totalb
KRT20_rs17881320_3Mbp_all.gemma.assoc.totalb
KRT20_rs17881320_3Mbp_noneczema.gemma.assoc.totalb
KRT33A_rs12951971_3Mbp_noneczema.gemma.assoc.totalb
KRT33A_rs17881320_3Mbp_noneczema.gemma.assoc.totalb
KRTAP1-4_rs12951971_3Mbp_all.gemma.assoc.totalb
KRTAP1-4_rs17881320_3Mbp_all.gemma.assoc.totalb
OVOL1_rs10791824_100kbp_all.gemma.assoc.totalb
OVOL1_rs10791824_100kbp_noneczema.gemma.assoc.totalb
OVOL1_rs10791824_1Mbp_all.gemma.assoc.totalb
OVOL1_rs10791824_1Mbp_eczema.gemma.assoc.totalb
OVOL1_rs10791824_1Mbp_noneczema.gemma.assoc.totalb
OVOL1_rs10791824_250kbp_all.gemma.assoc.totalb
OVOL1_rs10791824_250kbp_noneczema.gemma.assoc.totalb
OVOL1_rs10791824_3Mbp_all.gemma.assoc.totalb
OVOL1_rs10791824_3Mbp_eczema.gemma.assoc.totalb
OVOL1_rs10791824_3Mbp_noneczema.gemma.assoc.totalb
PGS1_rs11657987_1Mbp_all.gemma.assoc.totalb
PGS1_rs11657987_1Mbp_noneczema.gemma.assoc.totalb
PGS1_rs11657987_3Mbp_all.gemma.assoc.totalb
PGS1_rs11657987_3Mbp_noneczema.gemma.assoc.totalb
PPP2R3C_rs2038255_1Mbp_noneczema.gemma.assoc.totalb
PPP2R3C_rs2038255_3Mbp_noneczema.gemma.assoc.totalb
PPP2R3C_rs2143950_1Mbp_noneczema.gemma.assoc.totalb
PPP2R3C_rs2143950_3Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs12295535_1Mbp_all.gemma.assoc.totalb
PRR5L_rs12295535_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs12295535_3Mbp_all.gemma.assoc.totalb
PRR5L_rs12295535_3Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_1Mbp_all.gemma.assoc.totalb
PRR5L_rs2218565_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2218565_3Mbp_all.gemma.assoc.totalb
PRR5L_rs2218565_3Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_1Mbp_all.gemma.assoc.totalb
PRR5L_rs2592555_1Mbp_noneczema.gemma.assoc.totalb
PRR5L_rs2592555_3Mbp_all.gemma.assoc.totalb
PRR5L_rs2592555_3Mbp_noneczema.gemma.assoc.totalb
PRRC2A_rs145809981_250kbp_all.gemma.assoc.totalb
PRRC2A_rs41293864_250kbp_all.gemma.assoc.totalb
PSMA6_rs2038255_1Mbp_all.gemma.assoc.totalb
PSMA6_rs2038255_1Mbp_noneczema.gemma.assoc.totalb
PSMA6_rs2038255_3Mbp_all.gemma.assoc.totalb
PSMA6_rs2038255_3Mbp_noneczema.gemma.assoc.totalb
PSMA6_rs2143950_1Mbp_all.gemma.assoc.totalb
PSMA6_rs2143950_1Mbp_noneczema.gemma.assoc.totalb
PSMA6_rs2143950_3Mbp_all.gemma.assoc.totalb
PSMA6_rs2143950_3Mbp_noneczema.gemma.assoc.totalb
STAT6_rs1799986_1Mbp_all.gemma.assoc.totalb
STMN2_rs6473227_3Mbp_all.gemma.assoc.totalb
STMN2_rs6473227_3Mbp_noneczema.gemma.assoc.totalb
UGT3A1_rs10214237_1Mbp_all.gemma.assoc.totalb
UGT3A1_rs10214237_3Mbp_all.gemma.assoc.totalb


#Prepare Gassocplot2 input - separately for noneczema and all.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in *${my_snps}*1Mbp_*noneczema*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $2}')
snp_chrom=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $1}')
snp_pos=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_1Mbp_twinsuk_noneczema.gassocplot2
done
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in *${my_snps}*3Mbp_*noneczema*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $2}')
snp_chrom=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $1}')
snp_pos=$(grep -w $my_snp ${my_file%.totalb}.txt | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_3Mbp_twinsuk_noneczema.gassocplot2
done
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in *${my_snps}*3Mbp_*noneczema*totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
echo $snp $gene	$pval >>${snp}.3Mbp_noneczema.PPH4
done
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in *${my_snps}*1Mbp_*noneczema*totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
echo $snp $gene	$pval >>${snp}.1Mbp_noneczema.PPH4
done
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_3Mbp_twinsuk_noneczema.gassocplot2 | tail -1 > ${my_snps}_3Mbp_twinsuk_noneczema.top
	python $coloc_scripts/generate_gassoc_input2.py \
	${my_snps}_3Mbp_twinsuk_noneczema.top \
	${my_snps}.3Mbp_noneczema.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_3Mbp_twinsuk_skin_noneczema2
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_1Mbp_twinsuk_noneczema.gassocplot2 | tail -1 > ${my_snps}_1Mbp_twinsuk_noneczema.top
	python $coloc_scripts/generate_gassoc_input2.py \
	${my_snps}_1Mbp_twinsuk_noneczema.top \
	${my_snps}.1Mbp_noneczema.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_1Mbp_twinsuk_skin_noneczema2
done

eqtlgen=$HOME/analysis/colocalization/coloc/eqtlgen/euro/gene
sun=$HOME/analysis/colocalization/coloc/sun_pqtl

cd $analysis/gemma/plots
#Plot Sun2018 and eqtlgen, along with TwinsUK in one figure, if possible.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_3Mbp_sun2 ]
then
    Rscript $scripts/plot_stack_gassoc2_four.R \
    $analysis/gemma/skin/output/${my_snps}_3Mbp_twinsuk_skin_noneczema2 \
    $analysis/gemma/lcl/output/${my_snps}_3Mbp_twinsuk_lcl_noneczema2 \
    $eqtlgen/${my_snps}_3Mbp_eqtlgen2 \
    $sun/${my_snps}_3Mbp_sun2 \
	"eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen"  "pQTL blood - Sun (2018)"   ${my_snps}_3Mbp_stack_gassocplot4.pdf 
else
    Rscript $scripts/plot_stack_gassoc2_three.R \
    $analysis/gemma/skin/output/${my_snps}_3Mbp_twinsuk_skin_noneczema2 \
    $analysis/gemma/lcl/output/${my_snps}_3Mbp_twinsuk_lcl_noneczema2 \
    $eqtlgen/${my_snps}_3Mbp_eqtlgen2 \
    "eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen" \
    ${my_snps}_3Mbp_stack_gassocplot3.pdf
fi
done

#New edition of figures with gassocplot - this time, dot in the middle of gene, for all genes and labelled with the gene name.
#Missing datapoint not plotted.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	python $coloc_scripts/generate_gassoc_input3.py \
	${my_snps}_3Mbp_twinsuk_noneczema.top \
	${my_snps}.3Mbp_noneczema.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_3Mbp_twinsuk_skin_noneczema3
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	python $coloc_scripts/generate_gassoc_input3.py \
	${my_snps}_1Mbp_twinsuk_noneczema.top \
	${my_snps}.1Mbp_noneczema.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_1Mbp_twinsuk_skin_noneczema3
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_3Mbp_sun3 ]
then
    Rscript $scripts/plot_stack_gassoc2_four.R \
    $analysis/gemma/skin/output/${my_snps}_3Mbp_twinsuk_skin_noneczema3 \
    $analysis/gemma/lcl/output/${my_snps}_3Mbp_twinsuk_lcl_noneczema3 \
    $eqtlgen/${my_snps}_3Mbp_eqtlgen3 \
    $sun/${my_snps}_3Mbp_sun3 \
	"eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen"  "pQTL blood - Sun (2018)"   ${my_snps}_3Mbp_stack_gassocplot_3.pdf 
else
    Rscript $scripts/plot_stack_gassoc2_three.R \
    $analysis/gemma/skin/output/${my_snps}_3Mbp_twinsuk_skin_noneczema3 \
    $analysis/gemma/lcl/output/${my_snps}_3Mbp_twinsuk_lcl_noneczema3 \
    $eqtlgen/${my_snps}_3Mbp_eqtlgen3 \
    "eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen" \
    ${my_snps}_3Mbp_stack_gassocplot_3.pdf
fi
done

#Now, for 1Mbp interval.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_1Mbp_sun3 ]
then
    Rscript $scripts/plot_stack_gassoc2_four.R \
    $analysis/gemma/skin/output/${my_snps}_1Mbp_twinsuk_skin_noneczema3 \
    $analysis/gemma/lcl/output/${my_snps}_1Mbp_twinsuk_lcl_noneczema3 \
    $eqtlgen/${my_snps}_1Mbp_eqtlgen3 \
    $sun/${my_snps}_1Mbp_sun3 \
	"eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen"  "pQTL blood - Sun (2018)"   ${my_snps}_1Mbp_stack_gassocplot_3.pdf 
else
    Rscript $scripts/plot_stack_gassoc2_three.R \
    $analysis/gemma/skin/output/${my_snps}_1Mbp_twinsuk_skin_noneczema3 \
    $analysis/gemma/lcl/output/${my_snps}_1Mbp_twinsuk_lcl_noneczema3 \
    $eqtlgen/${my_snps}_1Mbp_eqtlgen3 \
    "eQTL skin - TwinsUK" "eQTL LCL - TwinsUK" "eQTL blood - eQTLGen" \
    ${my_snps}_1Mbp_stack_gassocplot_3.pdf
fi
done

