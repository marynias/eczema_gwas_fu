#!/bin/bash
###Replication of analysis by Mehsa Babaei and Ashley Budu-Aggrey. 
###This is analysis looking only at CD207 locus.
###Note: /projects/MRC-IEU is invisible from worker nodes.
#to identify co-localised GWAS and eQTL signal for locus on chr 2.
set -o nounset -o pipeail -o errexit
#1) Get list of SNPs for region 
#(eg, CD207 at chr2p13.3= chr2:71025000-71150000)
scripts=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/analysis/mehsa
utils=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/utils
analysis_dir=/panfs/panasas01/sscm/qh18484/analysis/mehsa
gxp_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk
var_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk

cd $analysis_dir

#your list of SNPs/potential eQTLs
CD207_snps="SNP_list.txt"
#genotype files (.gen and .sample files)
genfile=$var_data/gen/data.chr02.gen
sample_file=$analysis_dir/examples/chr2_Eurobats_Public_alt_pheno.sample
outfile=$analysis_dir/chr2_Eurobats_Public_subset_samples
exp_sample_id=$analysis_dir/examples/Sample_ID.txt
#number of genotyped individuals with expression data
number_of_samples=672
#bimbam file to be created from .gen file. Will need this to run gemma
mean_geno=$analysis_dir/chr2_genotype.bimbam
#modified sample file with expression data and age variables
new_sample_file=$analysis_dir/chr2_CD207_Eurobats_Public_exp_pheno.sample
# phenotype file with cases control status and expression data
pheno_file=$analysis_dir/chr2_pheno.bimbam
# covariate file with intercept column (1) and sex
covar_file=$analysis_dir/covar.txt
# file with SNP name (rs ID's), bp and chromosome

snp_file=$analysis_dir/examples/SNP_annotation_file.txt
#output file for eQTL analysis
eQTL_out=$analysis_dir/CD207_eQTL_analysis_test
#relationship matirx
relatedness_matrix=$analysis_dir/examples/relatedness_all_geno.cXX.txt
number_of_samples=672
merged_geno=$analysis_dir/All_chr_genotype.bimbam
temp_pheno=$analysis_dir/temp_pheno.txt

awk '$3 >= 71025000 && $3 <= 71150000 {print $2}' $genfile > $CD207_snps

#extract transcript data for CD207 (transcript id = ENSG00000116031, derived from GENCODE file)
sed -n -e '1p' -e '/ENSG00000116031/p' $gxp_data/skin/data.rpkm > CD207_transcript_data.txt

#Get the output from the script below - in its current form, it is not possible to run it without further data.
#Rscript $analysis_dir/scripts/get_expression_sample_id.R


#use sample id's to subset .gen and .sample imputation files
gtool -S --g $genfile --s $sample_file --og $outfile.gen --os $outfile.sample --sample_id $exp_sample_id --inclusion ${CD207_snps}

#Number of input samples: 807
#Samples...
#Number of output samples: 672
#Snp_id...
#Gen...
#Number of input SNPs: 534313
#Number of output SNPs: 448

#convert imputation .gen file to BIMBAM format (mean_genotype file):
cat $outfile.gen | awk -v s=${number_of_samples} \
'{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) \
printf "," $(i*3+3)*2+$(i*3+4); \
printf "\n" }' > ${mean_geno}

#create phenotype and covariate file: first merge the expression data
#This step below has not been done as requires inputs not present. 
#Just downloading its output and carrying on.
#Rscript $analysis_dir/scripts/merge_exp_data.R

#create phenotype file - extract the phenotype column from the sample file
#first column = case/control status (1 = case, 0 = controls)
#second column = expression data for CD207

tail -n +3 ${new_sample_file} | awk '{print $6, $7}' > ${pheno_file}

#create covariate file
#first column = intercept ("1"), second column = age

tail -n +3 ${new_sample_file} | awk '{print $9, $8}' > ${covar_file}

#create snp annotation file in R. Again, not using the script but just 
#taking the output - SNP_annotation_file.txt
#The analysis focused on 522 credible SNP set around CD207-related locus.
#Rscript $analysis_dir/scripts/GEMMA_snp_file.R

#####create relatedness matrix using all genotype data 
#will need to change sample_id's and no. of samples

#first subset individuals with exp data. This loop does not work - gen file is somehow not found, even though it works for script below.
##Note: the reason for this behaviour was because the genfile was in the MRC-IEU project folder, which is invisible from the nodes. The file path has been corrected since then. 
#for i in {01..22}
#do
#qsub -v InGen=$var_data/gen/data.chr${i}.gen,InSample=$var_data/sample/data.chr${i}.sample,OutGen=$analysis_dir/data_subset.chr${i}.gen,OutSample=$analysis_dir/data_subset.chr${i}.sample,Subset=$analysis_dir/examples/Sample_ID.txt,LogFile=$analysis_dir/subset_all_chr${i}.log $utils/sub_gtool_subset.sh
#done

bash $scripts/gtool_subset.sh

#convert imputation .gen file to BIMBAM format (mean_genotype file).
for i in {01..22}
do
qsub -v number_of_samples=672,gen_file=$analysis_dir/data_subset.chr${i}.gen,bb_file=$analysis_dir/chr${i}_gen.bimbam $utils/sub_gen2bimbam.sh
done

#merge into one file
cat $analysis_dir/chr*_gen.bimbam > $merged_geno

#create fake pheno_file where there is no missing phenotype data (else those samples will not be included when computing the relatedness matrix)
for i in $(seq 1 $number_of_samples); do echo 1; done > ${temp_pheno}

#create relatedness matrix. Warning! Gemma puts all the results into newly created output directory,
qsub $scripts/gemma.sh

#run eQTL analysis 
# "-n 2" will tell it to look at the second column of the pheno file that has the expression data
# use "-n 1 2" to look at both expression and case-control status
gemma -g $mean_geno -p $pheno_file -n 2 -a $snp_file -c $covar_file -lmm 4 -o CD207_eQTL_analysis_test -maf 0 -k $relatedness_matrix

#Sort by the P-value, Wald-adjusted
cat CD207_eQTL_analysis_test.assoc.txt | awk 'NR<2{print $0;next}{print $0| "sort -k 12 -g"}' >CD207_eQTL_analysis_test.assoc_sorted.txt

#Run colocalization analysis
Rscript --vanilla $scripts/run_coloc.R