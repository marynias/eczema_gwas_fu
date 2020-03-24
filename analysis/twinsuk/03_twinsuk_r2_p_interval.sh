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

cd $analysis/gemma/r2_interval
#Subset initial 3Mbp GWAS intervals to the intervals defined in the r2-based interval analysis using 1K European data.
#TwinsUK
for a in $analysis/rel/*_3Mbp.gwas
do
my_file=$(basename $a)
my_rsid=$(echo $my_file | cut -d"_" -f1)
my_start=$(grep -w $my_rsid $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed | cut -f2)
my_end=$(grep -w $my_rsid $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed | cut -f3)
python $utils/subset_interval.py --tab $a --header_tab Y --pos_tab 3 --int_start $my_start --int_end $my_end --output ${my_rsid}_r2_0.2_1k.gwas
done

#Run coloc - skin.
for a in $analysis/gemma/skin/output/*3Mbp_all.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 672 ${my_output}_skin_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

for a in $analysis/gemma/skin/output/*3Mbp_eczema.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 86 ${my_output}_skin_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

for a in $analysis/gemma/skin/output/*3Mbp_noneczema.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 530 ${my_output}_skin_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

#Run coloc - lcl
for a in $analysis/gemma/lcl/output/*3Mbp_all.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 764 ${my_output}_lcl_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

for a in $analysis/gemma/lcl/output/*3Mbp_eczema.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 92 ${my_output}_lcl_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

for a in $analysis/gemma/lcl/output/*3Mbp_noneczema.gemma.assoc.txt
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f2)
my_output=$(echo $my_file | cut -d"_" -f1,2,4 | cut -d "." -f1)
Rscript $scripts/twinsuk_gemma_coloc_input2.R 610 ${my_output}_lcl_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a 
done

#Check the colocalisation results.
function pvals {
my_tissue=$1
my_pop=$2
for a in *${my_pop}_${my_tissue}_r2_0.2_1k.totalb
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a >>${my_pop}_${my_tissue}.sig
done
}

for a in eczema noneczema all
do
pvals skin $a
pvals lcl $a
done

#Generate summary table with significant results.
function summary_table {
my_tissue=$1
my_pop=$2
while read line
do
set -- $line
my_file=$1
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $my_file | cut -d"_" -f1)
interval=$(echo $my_file | cut -d"_" -f2)
target_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${my_file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $my_file)
pval_p=$(awk 'NR == 7 {print $2}' ${my_file%b}p)
echo -e $snp'\t'$interval'\t'$gene'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_pop}_${my_tissue}.coloc_results
done < ${my_pop}_${my_tissue}.sig
}

for a in eczema noneczema all
do
summary_table skin $a
summary_table lcl $a
done

#Prepare PPH4 summary file.
function prepare_pph4
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in *_${my_snps}_noneczema_${tissue}_r2_0.2_1k.totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f1)
echo $my_snps $gene	$pval >>${my_snps}_twinsuk_noneczema_${tissue}_r2_0.2_1k.PPH4
done
done
}

for a in skin lcl
do
prepare_pph4 $a
done


function summary_table_all_results {
my_tissue=$1
echo -e "Gene\tLead_SNP_interval\tTissue\tTOTAL.PP.H4_pval\tTOTAL.PP.H4_beta\tTop_SNP\tTop_SNP.PP.H4" >>${my_tissue}.all.coloc_results
for file in *_all_${my_tissue}*.totalp
do
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
tissue=$(echo $my_file | cut -d"_" -f4- | sed 's/_r2_0.2_1k.totalp//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_b=$(awk 'NR == 7 {print $2}' ${file%p}b)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.all.coloc_results
done
}

for a in skin lcl
do
summary_table_all_results $a
done
