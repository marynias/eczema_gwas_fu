#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/figures
data_manipulation=$HOME/analysis/annotation/data_manipulation
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
twinsuk_scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
twinsuk=$HOME/analysis/twinsuk/gemma/r2_interval
gtex=$HOME/analysis/colocalization/coloc/gtex
eqtlgen=$HOME/analysis/colocalization/coloc/eqtlgen/euro/gene/r2_interval
sun=$HOME/analysis/colocalization/coloc/sun_pqtl/r2_interval
blueprint=$HOME/working/data/Datasets/eQTL/Chen2016/colocalisation
cd $analysis

#TwinsUK figures
function twinsuk_gassocpot_input
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
	python $scripts/generate_gassoc_input4.py \
	$twinsuk/${my_snps}_twinsuk_noneczema_${tissue}_r2_0.2_1k.PPH4 \
	$data_manipulation/gencode_ensembl.v19.annotation.middle \
	$chrom \
	${my_snps}_twinsuk_noneczema_${tissue}_r2_0.2_1k.gassocplot
done
}

for a in skin lcl
do
twinsuk_gassocpot_input $a
done

#Plot TwinsUK
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
    Rscript $scripts/plot_stack_gassoc2.R \
    ${my_snps}_twinsuk_noneczema_skin_r2_0.2_1k.gassocplot \
    ${my_snps}_twinsuk_noneczema_lcl_r2_0.2_1k.gassocplot \
	"eQTL skin - TwinsUK" "eQTL LCL - TwinsUK"  ${my_snps}_twinsuk_r2_0.2_1k.pdf 
done

#GTEx figures
#TwinsUK figures
function gtex_gassocpot_input
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
	python $scripts/generate_gassoc_input4.py \
	$gtex/${my_snps}_gtex_${tissue}_r2_0.2_1k.PPH4 \
	$data_manipulation/gencode_ensembl.v19.annotation.middle \
	$chrom \
	${my_snps}_gtex_${tissue}_r2_0.2_1k.gassocplot
done
}

for a in Skin_Not_Sun_Exposed_Suprapubic \
Spleen \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes
do
gtex_gassocpot_input $a
done

#Blueprint figures
function blueprint_gassocpot_input
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
	python $scripts/generate_gassoc_input4.py \
	$blueprint/${my_snps}_blueprint_${tissue}_r2_0.2_1k.PPH4 \
	$data_manipulation/gencode_ensembl.v19.annotation.middle \
	$chrom \
	${my_snps}_blueprint_${tissue}_r2_0.2_1k.gassocplot
done
}

for a in tcel mono neut
do
blueprint_gassocpot_input $a
done


#Plot GTEx
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
    Rscript $scripts/plot_stack_gassoc2_six.R \
    ${my_snps}_gtex_Skin_Not_Sun_Exposed_Suprapubic_r2_0.2_1k.gassocplot \
 	${my_snps}_gtex_Skin_Sun_Exposed_Lower_leg_r2_0.2_1k.gassocplot \
    ${my_snps}_gtex_Cells_Transformed_fibroblasts_r2_0.2_1k.gassocplot \
    ${my_snps}_gtex_Whole_Blood_r2_0.2_1k.gassocplot \
    ${my_snps}_gtex_Cells_EBV-transformed_lymphocytes_r2_0.2_1k.gassocplot \
    ${my_snps}_gtex_Spleen_r2_0.2_1k.gassocplot \
	"eQTL Skin_Not_Sun_Exposed_Suprapubic - GTEx" \
	"eQTL Skin_Sun_Exposed_Lower_leg - GTEx" \
	"eQTL Cells_Transformed_fibroblasts - GTEx" \
	"eQTL Whole_Blood - GTEx" \
	"eQTL Cells_EBV-transformed_lymphocytes - GTEx" \
	"eQTL Spleen - GTEx" \
	${my_snps}_gtex_r2_0.2_1k.pdf 
done

#Prepare gassocplot input for Sun and Eqtlgen.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
	python $scripts/generate_gassoc_input4.py \
	$eqtlgen/${my_snps}_eqtlgen_r2_0.2_1k.PPH4 \
	$data_manipulation/gencode_ensembl.v19.annotation.middle \
	$chrom \
	${my_snps}_eqtlgen_r2_0.2_1k.gassocplot
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
	python $scripts/generate_gassoc_input4.py \
	$sun/${my_snps}_sun_r2_0.2_1k.PPH4 \
	$data_manipulation/gencode_ensembl.v19.annotation.middle \
	$chrom \
	${my_snps}_sun_r2_0.2_1k.gassocplot
done


#Plot TwinsUK, Sun and Eqtlgen together.
#Plot TwinsUK
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
    Rscript $twinsuk_scripts/plot_stack_gassoc2_four.R \
    ${my_snps}_twinsuk_noneczema_skin_r2_0.2_1k.gassocplot \
    ${my_snps}_twinsuk_noneczema_lcl_r2_0.2_1k.gassocplot \
    ${my_snps}_eqtlgen_r2_0.2_1k.gassocplot \
    ${my_snps}_sun_r2_0.2_1k.gassocplot \
	"eQTL skin - TwinsUK" \
	"eQTL LCL - TwinsUK" \
	"eQTL blood - eQTLGen" \
	"pQTL blood - Sun (2018)" \
	${my_snps}_twinsuk_eqtlgen_sun_r2_0.2_1k.pdf 
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_1Mbp_sun3 ]
then
    Rscript $twinsuk_scripts/plot_stack_gassoc2_four.R \
    ${my_snps}_twinsuk_noneczema_skin_r2_0.2_1k.gassocplot \
    ${my_snps}_twinsuk_noneczema_lcl_r2_0.2_1k.gassocplot \
    ${my_snps}_eqtlgen_r2_0.2_1k.gassocplot \
    ${my_snps}_sun_r2_0.2_1k.gassocplot \
	"eQTL skin - TwinsUK" \
	"eQTL LCL - TwinsUK" \
	"eQTL blood - eQTLGen" \
	"pQTL blood - Sun (2018)" \
	${my_snps}_twinsuk_eqtlgen_sun_r2_0.2_1k.pdf 
else
    Rscript $twinsuk_scripts/plot_stack_gassoc2_three.R \
    ${my_snps}_twinsuk_noneczema_skin_r2_0.2_1k.gassocplot \
    ${my_snps}_twinsuk_noneczema_lcl_r2_0.2_1k.gassocplot \
    ${my_snps}_eqtlgen_r2_0.2_1k.gassocplot \
	"eQTL skin - TwinsUK" \
	"eQTL LCL - TwinsUK" \
	"eQTL blood - eQTLGen" \
	${my_snps}_twinsuk_eqtlgen_r2_0.2_1k.pdf 
fi
done

#Plot Blueprint
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
    Rscript $twinsuk_scripts/plot_stack_gassoc2_three.R \
    ${my_snps}_blueprint_tcel_r2_0.2_1k.gassocplot \
    ${my_snps}_blueprint_mono_r2_0.2_1k.gassocplot \
    ${my_snps}_blueprint_neut_r2_0.2_1k.gassocplot \
	"eQTL T cells - Blueprint" \
	"eQTL monocytes - Blueprint" \
	"eQTL neutrophils - Blueprint" \
	${my_snps}_blueprint_r2_0.2_1k.pdf 
done

#Separate the results into a folder for each SNP.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12295535 rs11923593 rs2143950 rs17881320
do
	mkdir $my_snps
	cp ${my_snps}*pdf ./${my_snps}/
done