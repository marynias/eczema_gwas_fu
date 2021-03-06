#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/cedar
scripts=$HOME/bin/eczema_gwas_fu/analysis/cedar
coloc_scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
raw_data=$HOME/working/data/Datasets/eQTL/Momozawa2018
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
data_cleanup=$HOME/analysis/annotation/data_manipulation
twinsuk_scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk

cd $analysis
cd $raw_data

#Gene expression values corrected for sex, age, smoking, batch and number of Principal Components maximizing the number of cis-eQTLs in PLINK format 
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/read.me
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/Cedar_Pheno.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/read.me
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/CD15_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/TR_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/CD4_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/RE_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/PLA_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/IL_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/CD8_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/CD14_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GE/GE_Corr_4PCs_Cis/CD19_GE_Corrected4_Covars_PCs.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/read.me
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/SNPsRsIdMapping.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/CEDAR_Genotypes.bim.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/CEDAR_Genotypes.fam.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/CEDAR_Genotypes.bed.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/SNPsLoc_all.txt.gz
wget http://139.165.108.18/srv/genmol/permanent/1be6993fe41c12a051c9244d67c91da2be49e5dd26a6cd79f442bc006971e2ef/CEDAR_GENO/CEDAR_Imputed.dos.gz

#Check if the Plink files contain imputed genotypes. Check successful.
plink --bfile CEDAR_Genotypes --missing --recode oxford

cd $analysis

##Gene expression used - corrected for batch, sex, age, smoking and for principal components.
##Convert from Dosage MatrixEQTL format to Bimbam, via Geno.
#Convert Dosage MatrixEQTL format to Oxford Geno.
python $scripts/cedar2oxford.py $raw_data/SNPsRsIdMapping.txt $raw_data/CEDAR_Imputed.dos
#Copy the sample file generated by Plink into a file for each chromosome.
for a in {1..22}
do
cat $raw_data/plink.sample >cedar.chr${a}.sample 
done

#Generate a BIM file of the markers in the analysis.
for a in *.gen
do
cat $a | awk -v OFS="\t" '{print $1, $2, 0, $3, $4, $5}' >>cedar.bim
done

#Generate a CEDAR version of the GWAS results.
python $utils/update_rsid.py --bim cedar.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >results.euro.pval.cedar

#Subset GWAS results to the intervals defined in the r2-based interval analysis using 1K European data.
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d" " -f4)
python $utils/subset_interval_chr.py --tab results.euro.pval.cedar --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_r2_0.2_1k.gwas
done < $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed

#Add MAF column
for a in *gwas
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 

#Subset genotypes to SNPs of interest.
#Create BED files with intervals to select for each index SNP.
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-1500000, $3+1500000, $1}' >paternoster_2015_index_snps_sorted_3Mbp.bed
#Convert bim to bed.
cat cedar.bim | awk -v OFS="\t" '{print $1, $4, $4, $2}' >cedar.bed
while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_3Mbp.bed
done < paternoster_2015_index_snps_sorted_3Mbp.bed

#Identify SNPs overlapping the interval
for a in rs*.bed
do
intersectBed -wa -a cedar.bed -b $a>${a%.bed}_snps.bed
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

#Substitute Illumina IDs with Ensembl gene ids in the expression file.
#First, substitute gene names with Ensembl gene ids.

#Going to carry out analysis on the 3Mbp interval for now. 
#extract transcript data.
while read line
do
set -- $line
gene_id=$1
grep -w "$gene_id" $skin_data/data.rpkm >>skin/${b%.gene_ids}.rpkm
done < $b

for a in $raw_data/CD14_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD15_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD19_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD4_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD8_GE_Corrected4_Covars_PCs.txt \
$raw_data/IL_GE_Corrected4_Covars_PCs.txt \
$raw_data/PLA_GE_Corrected4_Covars_PCs.txt \
$raw_data/RE_GE_Corrected4_Covars_PCs.txt \
$raw_data/TR_GE_Corrected4_Covars_PCs.txt
do
python $scripts/probe_to_eng.py $raw_data/Probes_good_reanno_31137__14c2355106e13d3ec50d_.txt $a $data_cleanup/hugo_synonyms_ids2_filtered.ensembl.upper > ${a%.txt}.eng.txt
#Transpose the table with expression data.
awk '{for (i=1; i<=NF; i++) a[i,NR]=$i; max=(max<NF?NF:max)} END {for (i=1; i<=max; i++) {for (j=1; j<=NR; j++) printf "%s%s", a[i,j], (j==NR?RS:FS) }}' ${a%.txt}.eng.txt >${a%.txt}.eng.transposed.txt
done

for a in $raw_data/CD14_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD15_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD19_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD4_GE_Corrected4_Covars_PCs.txt \
$raw_data/CD8_GE_Corrected4_Covars_PCs.txt \
$raw_data/IL_GE_Corrected4_Covars_PCs.txt \
$raw_data/PLA_GE_Corrected4_Covars_PCs.txt \
$raw_data/RE_GE_Corrected4_Covars_PCs.txt \
$raw_data/TR_GE_Corrected4_Covars_PCs.txt
do
my_input=${a%.txt}.eng.transposed.txt
my_base=$(basename "$my_input")
sample=$(echo $my_base | cut -d"." -f1)
for b in rs*3Mbp_snps.gene_ids
do
while read line
do
set -- $line
gene_id=$1
grep -w "$gene_id" $my_input >>${b%_snps.gene_ids}_${sample}.expr
done < $b
done
done

#Fix header in each file.
for a in *.expr
do
my_sample=$(echo $a | cut -d"_" -f3- | sed 's/.expr//')
my_initial=${my_sample}.eng.transposed.txt
head -1 $raw_data/$my_initial | cat - $a >${a%.expr}_headers.expr
done

#####Subset genotype file.
#Create files with SNP ids in each interval.
for my_file in rs*_snps.bed
do
cut -f4 $my_file > ${my_file%.bed}.lst
done

#Create a list of individuals.
cut -d" " -f1 $raw_data/plink.sample | tail -n +3 >expression_individuals

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
expression_file=$raw_data/${sample}_GE_Corrected4_Covars_PCs.eng.transposed.txt
head -1 $expression_file | tr ' ' '\n' | tail -n +2 >${expression_file%.txt}.individuals
gtool -S --g cedar.chr${chrom}.gen --s cedar.chr${chrom}.sample --og ${snp}_3Mbp_${sample}_cedar.gen \
--os ${snp}_3Mbp_${sample}_cedar.sample --sample_id ${expression_file%.txt}.individuals --inclusion ${snp}_3Mbp_snps.lst
done
done


#Filter out SNPs with no rsid label.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
awk '$2 != "."' ${snp}_3Mbp_${sample}_cedar.gen >temp
mv temp ${snp}_3Mbp_${sample}_cedar.gen
done
done

#Convert gen files to bimbam.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
my_samples=$(wc -l $raw_data/${sample}_GE_Corrected4_Covars_PCs.eng.transposed.individuals | cut -d" " -f1)
qsub -v gen_file=${snp}_3Mbp_${sample}_cedar.gen,number_of_samples=$my_samples,bb_file=${snp}_3Mbp_${sample}_cedar.bimbam $utils/sub_gen2bimbam.sh
done
done

#Prepare SNP annotation file (SNP id, SNP position and SNP chrom)
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
python $twinsuk_scripts/gemma_snp_annotation.py cedar.bim ${snp}_3Mbp_${sample}_cedar.bimbam
done
done

#list of individuals in the gen file
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
cat ${snp}_3Mbp_${sample}_cedar.sample | cut -d" " -f2 | tail -n +3 >${snp}_3Mbp_${sample}_cedar.sample.lst
done
done

#Preparing relatedness matrix.
#Merge all the BIMBAM files into one matrix.
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
cat *_3Mbp_${sample}_cedar.bimbam > ${sample}_3Mbp_cedar.bimbam 
done

#create fake pheno_file where there is no missing phenotype data (else those samples will not be included when computing the relatedness matrix)
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
no_inds=$(wc -l rs7512552_3Mbp_${sample}_cedar.sample.lst | cut -d" " -f1)
for i in $(seq 1 $no_inds) ; do echo 1; done > ${sample}_3Mbp_cedar.pheno
done

for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
qsub -v my_bimbam=${sample}_3Mbp_cedar.bimbam,my_pheno=${sample}_3Mbp_cedar.pheno,my_output=${sample}_3Mbp_cedar.rel $scripts/sub_gemma_rel.sh
done

##Create gene expression files.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
Rscript --vanilla $scripts/cedar_sort_expression.R ensemble_ids_unique_ensembl_names ${snp}_3Mbp_${sample}_GE_Corrected4_Covars_PCs_headers.expr \
${snp}_3Mbp_${sample}_cedar.sample.lst
done
done

#Run Gemma.
for file in *gxp
do
gene=$(echo $file | cut -d"_" -f1)
snp=$(echo $file | cut -d"_" -f2)
sample=$(echo $file | cut -d"_" -f4 | cut -d"." -f1)
gemma -g ${snp}_3Mbp_${sample}_cedar.bimbam -p $file -a ${snp}_3Mbp_${sample}_cedar.snps -k output/${sample}_3Mbp_cedar.rel.cXX.txt \
-lmm 4 -o ${gene}_${snp}_3Mbp_${sample}.gemma
done

#Run coloc.
for file in output/*.gemma.assoc.txt
do
snp=$(echo $file | cut -d"_" -f2)
sample=$(echo $file | cut -d"_" -f4 | cut -d"." -f1)
no_inds=$(wc -l ${snp}_3Mbp_${sample}_cedar.sample.lst | cut -d" " -f1)
Rscript $scripts/cedar_coloc_input.R ${snp}_r2_0.2_1k.gwas $file $no_inds
done

#Check the colocalisation results.
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
for a in *_${sample}.gemma.assoc.totalb
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a >>cedar_${sample}_significant_totalb.txt
done
done

function summary_table {
my_tissue=$1
while read line
do
set -- $line
file=$1
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $file)
pval_p=$(awk 'NR == 7 {print $2}' ${file%b}p)
echo -e $gene'\t'$snp'\t'$my_tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>cedar_${my_tissue}.coloc_results.b
done < cedar_${my_tissue}_significant_totalb.txt
}

for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
summary_table $sample
done


function summary_table_all_results {
my_tissue=$1
echo -e "Gene\tLead_SNP_interval\tTissue\tTOTAL.PP.H4_pval\tTop_SNP\tTop_SNP.PP.H4" >>${my_tissue}.all.coloc_results
for file in *${my_tissue}*.totalp
do
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
tissue=$(echo $my_file | cut -d"_" -f4- | sed 's/.gemma.assoc.totalp//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.all.coloc_results
done
}

for a in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
summary_table_all_results $a
done
