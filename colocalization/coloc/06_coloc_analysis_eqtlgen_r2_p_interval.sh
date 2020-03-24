#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/eqtlgen/euro/gene/r2_interval
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
eqtlgen=/panfs/panasas01/sscm/qh18484/data/eqtl/eQTLgen
utils=$HOME/bin/eczema_gwas_fu/utils
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
eqtl_f=eQTLsFDR-ProbeLevel.txt.gz
gwas_f=results.euro.pval.eqtlgen.ea

cd $analysis

#Subset initial 3Mbp GWAS intervals to the intervals defined in the r2-based interval analysis using 1K European data.
#Eqtlgen
for a in ../rs*_3Mbp.gwas 
do
my_file=$(basename $a)
my_rsid=$(echo $my_file | cut -d"_" -f1)
my_start=$(grep -w $my_rsid $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed | cut -f2)
my_end=$(grep -w $my_rsid $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed | cut -f3)
python $utils/subset_interval.py --tab $a --header_tab Y --pos_tab 3 --int_start $my_start --int_end $my_end --output ${my_rsid}_r2_0.2_1k.gwas
done

#Run coloc analysis
for a in ../*3Mbp_*.eqtl
do
my_file=$(basename $a)
my_gwas=$(echo $my_file | cut -d"_" -f1)
my_output=$(echo $my_file | cut -d"_" -f1,3 | cut -d"." -f1)
echo ${my_gwas}_r2_0.2_1k.gwas
echo $a 
Rscript $scripts/coloc_eqtlgen3.R ${my_output}_r2_0.2_1k ${my_gwas}_r2_0.2_1k.gwas $a 
done 

#Check the colocalisation results.
for a in *.totalb
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>eqtlgen.sig
done

#Generate summary table with significant results.
while read line
do
set -- $line
my_file=$1
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
target_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${my_file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $my_file)
pval_p=$(awk 'NR == 7 {print $2}' ${my_file%b}p)
echo -e $gene'\t'$snp'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>eqtlgen.coloc_results
done < eqtlgen.sig

#Prepare PPH4 summary file.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*_r2_0.2_1k.totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f2 | sed 's/_r2_0.2_1k.totalb//')
echo $my_snps $gene	$pval >>${my_snps}_eqtlgen_r2_0.2_1k.PPH4
done
done

function summary_table_all_results {
echo -e "Gene\tLead_SNP_interval\tTOTAL.PP.H4_pval\tTOTAL.PP.H4_beta\tTop_SNP\tTop_SNP.PP.H4" >>eqtlgen.all.coloc_results
for file in *.totalp
do
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_b=$(awk 'NR == 7 {print $2}' ${file%p}b)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>eqtlgen.all.coloc_results
done
}

summary_table_all_results 


