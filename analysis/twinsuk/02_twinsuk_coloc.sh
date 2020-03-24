#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/twinsuk
coloc=$HOME/analysis/twinsuk/rel
coloc_scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk
skin_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/skin
lcl_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/LCLs
var_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/gen
sample_data=/panfs/panasas01/sscm/qh18484/data/eqtl/twinsuk/sample
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils

#Colocalization analysis
cd $coloc

#Extract gene names from Ensembl identifiers.
cd $analysis/skin
Rscript --vanilla $scripts/extract_list_ids_ensembl.R ensemble_ids_unique

cd $analysis/lcl
Rscript --vanilla $scripts/extract_list_ids_ensembl.R ensemble_ids_unique

#Extract also transcript end, transcript start and the middle point of gene.


#Modify the REL matrixeqtl output to calculate SE, Z-score and change Ensembl Gene identifiers to gene names.
function calculate {
my_dir=$1
cd $my_dir
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp 
do
Rscript --vanilla $scripts/calculate_matrixeqtl.R ${snp}_${interval}_all_rel.matrixeqtl.matrixeqtl ensemble_ids_unique_ensembl_names
Rscript --vanilla $scripts/calculate_matrixeqtl.R ${snp}_${interval}_eczema_rel.matrixeqtl.matrixeqtl ensemble_ids_unique_ensembl_names
Rscript --vanilla $scripts/calculate_matrixeqtl.R ${snp}_${interval}_noneczema_rel.matrixeqtl.matrixeqtl ensemble_ids_unique_ensembl_names
done
done
}

calculate $analysis/skin
calculate $analysis/lcl

gwas_f=$gwas/results.euro.pval.twinsuk.harmonized_beta
#Export GWAS results in a given interval.
for my_int in 100 250 1000 3000
do
python $coloc_scripts/generate_gwas_input.py --tab $gwas_f \
--proces $gwas/paternoster_2015_index_snps_sorted.txt --interval $my_int --chrom 2 --pos 3 --eaf 7 --ident 12 
done

#Rename GWAS files from kbp to Mbp
for a in *000kbp.gwas
do
mv $a ${a%000kbp.gwas}Mbp.gwas
done


#Prepare input from Twinsuk eQTL results.
function input_eqtl {
my_tissue=$1
cd $analysis/$my_tissue
for snp in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for interval in 1Mbp 3Mbp 100kbp 250kbp
do
python $scripts/create_eqtl_input_matrixeqtl.py ${snp}_${interval}_all_rel.matrixeqtl_final $analysis/twinsuk.bed $coloc/$my_tissue/${snp}_${interval}_all_rel.matrixeqtl.eqtl
python $scripts/create_eqtl_input_matrixeqtl.py ${snp}_${interval}_eczema_rel.matrixeqtl_final $analysis/twinsuk.bed $coloc/$my_tissue/${snp}_${interval}_eczema_rel.matrixeqtl.eqtl
python $scripts/create_eqtl_input_matrixeqtl.py ${snp}_${interval}_noneczema_rel.matrixeqtl_final $analysis/twinsuk.bed $coloc/$my_tissue/${snp}_${interval}_noneczema_rel.matrixeqtl.eqtl
done
done
}

input_eqtl skin
input_eqtl lcl

#Divide the resulting master file for gene-by-gene files.
function create_eqtl_gene {
my_tissue=$1
cd $coloc/$my_tissue
for infile in *matrixeqtl.eqtl 
do
qsub -v my_file=$infile $scripts/sub_divide_by_gene.sh
done
}

create_eqtl_gene skin
create_eqtl_gene lcl

#Convert from spaces to tabs and add header
function fix_eqtl {
my_tissue=$1
cd $coloc/$my_tissue
for my_input in *matrixeqtl_*.eqtl 
do
sed -i '1s/^/rsid\tchrom\tpos\tpval\tbeta\tse\tgene\n/' $my_input
awk -v OFS='\t' '{$1=$1;print}' $my_input >temp 
mv temp $my_input
done
}

fix_eqtl skin
fix_eqtl lcl

function coloc_all {
my_tissue=$1
my_samples=$2
cd $coloc/$my_tissue
for my_input in *all*.matrixeqtl_*.eqtl  
do
echo $my_input
my_gwas=$(echo $my_input |cut -d"_" -f1,2)
Rscript $scripts/twinsuk_coloc_input.R $coloc/${my_gwas}.gwas $my_input $my_samples
done
}

function coloc_eczema {
my_tissue=$1
my_samples=$2
cd $coloc/$my_tissue
for my_input in *_eczema_*.matrixeqtl_*.eqtl  
do
echo $my_input
my_gwas=$(echo $my_input |cut -d"_" -f1,2)
Rscript $scripts/twinsuk_coloc_input.R $coloc/${my_gwas}.gwas $my_input $my_samples
done
}

function coloc_noneczema {
my_tissue=$1
my_samples=$2
cd $coloc/$my_tissue
for my_input in *_noneczema_*.matrixeqtl_*.eqtl  
do
echo $my_input
my_gwas=$(echo $my_input | cut -d"_" -f1,2)
Rscript $scripts/twinsuk_coloc_input.R $coloc/${my_gwas}.gwas $my_input $my_samples
done
}

coloc_all skin 672
coloc_eczema skin 86
coloc_noneczema skin 530

coloc_all lcl 764
coloc_eczema lcl 92
coloc_noneczema lcl 610

#Identify the genes with overall significant PP.H4 calculated using p-vals.
function pvals {
my_tissue=$1
cd $coloc/$my_tissue
for a in *totalp
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done

for a in *totalb
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done
}

pvals skin
pvals lcl

#Skin significant results
rs12951971_1Mbp_all_rel.matrixeqtl_DNAJC7.totalb
rs17881320_1Mbp_all_rel.matrixeqtl_DNAJC7.totalb

#LCL significant results
rs10791824_100kbp_all_rel.matrixeqtl_OVOL1.totalp
rs10791824_250kbp_all_rel.matrixeqtl_OVOL1.totalp
rs12951971_3Mbp_all_rel.matrixeqtl_KRTAP9-1.totalp
rs17881320_3Mbp_all_rel.matrixeqtl_KRTAP9-1.totalp
rs2918307_100kbp_noneczema_rel.matrixeqtl_ACTL9.totalp
rs2918307_250kbp_noneczema_rel.matrixeqtl_ACTL9.totalp
rs10214237_100kbp_all_rel.matrixeqtl_IL7R.totalb
rs10791824_100kbp_all_rel.matrixeqtl_OVOL1.totalb
rs10791824_250kbp_all_rel.matrixeqtl_OVOL1.totalb
rs12188917_1Mbp_all_rel.matrixeqtl_RAD50.totalb
rs2918307_100kbp_noneczema_rel.matrixeqtl_ACTL9.totalb
rs2918307_250kbp_noneczema_rel.matrixeqtl_ACTL9.totalb
rs4705962_1Mbp_all_rel.matrixeqtl_RAD50.totalb

function create_gassocplot_input {
my_tissue=$1
my_cohort=$2
my_interval=$3
cd $coloc/$my_tissue
#Prepare PPH4 input for plotting with gassocplot -
for my_snps in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_file in ${my_snps}_${my_interval}_${my_cohort}_rel*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $1}')
snp_chrom=$(grep $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $2}')
snp_pos=$(grep $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval >>${my_snps}_${my_interval}_${my_cohort}_rel_twinsuk.gassocplot
done
done
}

for a in skin #lcl
do
for b in eczema noneczema all
do
create_gassocplot_input $a $b 1Mbp
create_gassocplot_input $a $b 3Mbp
done
done


#Divide the results into one directory for SNPs, each.
function create_dirs {
my_tissue=$1
for my_snps in rs112111458 rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 1Mbp 3Mbp 100kbp 250kbp
do
cd $coloc/$my_tissue
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
mkdir -p chr${chrom}_${my_snps}/${my_int}
cp -r ${my_snps}*${my_int}* chr${chrom}_${my_snps}/${my_int}
done
done
}

create_dirs skin
create_dirs lcl