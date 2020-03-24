#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/eqtlgen/euro/gene
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
eqtlgen=/panfs/panasas01/sscm/qh18484/data/eqtl/eQTLgen
utils=$HOME/bin/eczema_gwas_fu/utils
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
eqtl_f=eQTLsFDR-ProbeLevel.txt.gz
gwas_f=results.euro.pval.eqtlgen.ea

cd $analysis

#Gene by gene analysis of the eQTLgen dataset of genes within a given interval of focal SNP.

#Substitute the relevant ENST (transcript) ids with ENG (gene) ids.
#Annotate all the eczema target genes with ENsembl IDs
Rscript $scripts/match_gene_names2ensembl_ids.R $HOME/analysis/colocalization/coloc/sun_pqtl/all_targets_sun2018_3Mbp.uniq all_targets_sun2018_3Mbp.ensembl

#Substitute transcript IDS in eqtlgen files with gene names within 3MBbp of index SNP in eczema GWAS.
#Print a subset of eQTLgen dataset containing onlyla substituted lines.
python $scripts/eqtlgen_ensembl_ids.py all_targets_sun2018_3Mbp.ensembl $eqtlgen/eQTLsFDR-ProbeLevel.txt.gz QTLsFDR-ProbeLevel_genes_eczema3Mbp.txt

#Generate Euro GWAS file with rsids matching those from the eQTLgen study.
python $utils/update_rsid.py --bim $HOME/analysis/colocalization/coloc/eqtlgen/published/diagnostics/${eqtl_f%.txt.gz}_unique.bim \
--tab $gwas/results.euro.pval.tsv --head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$analysis/results.euro.pval.eqtlgen

#Harmonize the effect sizes so that they are relevant to the same allele in both files.
python $utils/harmonize_beta.py --tab $analysis/results.euro.pval.eqtlgen --ref $eqtlgen/$eqtl_f --header_tab Y --header_ref Y \
--rsid_tab 12 --rsid_ref 2 --effect_tab 4 --alt_tab 5 --effect_ref 10 --beta_tab 8 --zscore_tab 10 --out $analysis/$gwas_f

#Rename the Mbase intervals to Mbp, rather than Kbp.
for a in *000kbp.gwas
do
mv $a ${a%000kbp.gwas}Mbp.gwas
done

#Export eQTL results in a given interval.
#for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
#do
#for my_int in 10kbp 100kbp 250kbp 1Mbp 3Mbp
#do
#python $scripts/gene_generate_eqtl_input_pvals_zscore_by_rsid.py --tab $analysis/QTLsFDR-ProbeLevel_genes_eczema3Mbp.txt --gene 17 \
#--proc $gwas/paternoster_2015_index_snps_sorted.txt --interval $my_int --snp $my_snps \
#--genes $HOME/analysis/colocalization/coloc/sun_pqtl/${my_snps}_${my_int}.gene_names_abbrv --chrom 3 --pos 4 --pval 1 --ident 2 --zscore 11
#done
#done

for my_s in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_i in 10kbp 100kbp 250kbp 1Mbp 3Mbp
do
qsub -v my_snps=$my_s,my_int=$my_i $scripts/sub_generate_eqtl_input_eqtlgen.sh
done
done

#Divide the resulting master file for gene-by-gene files.
for my_file in *bp.eqtl
do
qsub -v a=$my_file $scripts/sub_process_eqtl_input_eqtlgen.sh
done

#Convert from spaces to tabs and add header
for my_input in *bp_*.eqtl 
do
sed -i '1s/^/rsid\tchrom\tpos\tpval\tzscore\tgene\n/' $my_input
awk -v OFS='\t' '{$1=$1;print}' $my_input >temp 
mv temp $my_input
done

#Run coloc analysis
for my_input in *bp_*.eqtl
do
my_gwas=$(echo $my_input | cut -d"_" -f1,2)
Rscript $scripts/coloc_eqtlgen2.R ${my_gwas}.gwas $my_input 
done

#Identify the genes with overall significant PP.H4 calculated using p-vals.
for a in *totalp
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a
done
#17 genes total.
rs10199605_3Mbp_ADAM17.totalp
rs10199605_3Mbp_IAH1.totalp
rs12295535_100kbp_PRR5L.totalp
rs12295535_10kbp_PRR5L.totalp
rs12295535_1Mbp_PRR5L.totalp
rs12295535_250kbp_PRR5L.totalp
rs12295535_3Mbp_PRR5L.totalp
rs12730935_100kbp_IL6R.totalp
rs12730935_10kbp_IL6R.totalp
rs12730935_1Mbp_IL6R.totalp
rs12730935_250kbp_IL6R.totalp
rs12730935_3Mbp_IL6R.totalp
rs12951971_250kbp_STAT5B.totalp
rs1799986_100kbp_LRP1.totalp
rs1799986_100kbp_STAT6.totalp
rs1799986_10kbp_LRP1.totalp
rs1799986_250kbp_LRP1.totalp
rs2038255_100kbp_FAM177A1.totalp
rs2038255_100kbp_PPP2R3C.totalp
rs2038255_10kbp_FAM177A1.totalp
rs2038255_10kbp_PPP2R3C.totalp
rs2038255_1Mbp_FAM177A1.totalp
rs2038255_1Mbp_PPP2R3C.totalp
rs2038255_250kbp_FAM177A1.totalp
rs2038255_250kbp_PPP2R3C.totalp
rs2038255_3Mbp_FAM177A1.totalp
rs2038255_3Mbp_PPP2R3C.totalp
rs2041733_10kbp_CLEC16A.totalp
rs2143950_100kbp_FAM177A1.totalp
rs2143950_100kbp_PPP2R3C.totalp
rs2143950_10kbp_FAM177A1.totalp
rs2143950_10kbp_PPP2R3C.totalp
rs2143950_1Mbp_FAM177A1.totalp
rs2143950_1Mbp_PPP2R3C.totalp
rs2143950_250kbp_FAM177A1.totalp
rs2143950_250kbp_PPP2R3C.totalp
rs2143950_3Mbp_FAM177A1.totalp
rs2143950_3Mbp_PPP2R3C.totalp
rs2212434_1Mbp_LRRC32.totalp
rs2212434_250kbp_LRRC32.totalp
rs2212434_3Mbp_LRRC32.totalp
rs2218565_100kbp_PRR5L.totalp
rs2218565_10kbp_PRR5L.totalp
rs2218565_1Mbp_PRR5L.totalp
rs2218565_250kbp_PRR5L.totalp
rs2218565_3Mbp_PRR5L.totalp
rs2592555_100kbp_PRR5L.totalp
rs2592555_10kbp_PRR5L.totalp
rs2592555_1Mbp_PRR5L.totalp
rs2592555_250kbp_PRR5L.totalp
rs2592555_3Mbp_PRR5L.totalp
rs3917265_1Mbp_IL18RAP.totalp
rs3917265_1Mbp_IL1RL1.totalp
rs3917265_3Mbp_IL18RAP.totalp
rs3917265_3Mbp_IL1RL1.totalp
rs4809219_1Mbp_LIME1.totalp
rs4809219_250kbp_LIME1.totalp
rs4809219_3Mbp_LIME1.totalp
rs6419573_100kbp_IL18RAP.totalp
rs6419573_1Mbp_IL18RAP.totalp
rs6419573_1Mbp_IL1RL1.totalp
rs6419573_250kbp_IL18RAP.totalp
rs6419573_250kbp_IL1RL1.totalp
rs6419573_3Mbp_IL18RAP.totalp
rs6419573_3Mbp_IL1RL1.totalp
rs6602364_100kbp_IL2RA.totalp
rs6602364_1Mbp_IL2RA.totalp
rs6602364_250kbp_IL2RA.totalp
rs6602364_3Mbp_IL2RA.totalp
rs7127307_1Mbp_ETS1.totalp
rs7127307_3Mbp_ETS1.totalp
rs7146581_100kbp_TRAF3.totalp
rs7146581_10kbp_TRAF3.totalp
rs7146581_1Mbp_TRAF3.totalp
rs7146581_250kbp_TRAF3.totalp
rs7146581_3Mbp_TRAF3.totalp

#Identify the genes with overall significant PP.H4 calculated using betas and variance of beta.
#15 genes total.
for a in *totalb
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a
done
rs10199605_3Mbp_ADAM17.totalb
rs10199605_3Mbp_IAH1.totalb
rs12295535_100kbp_PRR5L.totalb
rs12295535_10kbp_PRR5L.totalb
rs12295535_1Mbp_PRR5L.totalb
rs12295535_250kbp_PRR5L.totalb
rs12295535_3Mbp_PRR5L.totalb
rs12730935_100kbp_IL6R.totalb
rs12730935_10kbp_IL6R.totalb
rs12730935_1Mbp_IL6R.totalb
rs12730935_250kbp_IL6R.totalb
rs12730935_3Mbp_IL6R.totalb
rs12951971_250kbp_STAT5B.totalb
rs1799986_100kbp_LRP1.totalb
rs1799986_10kbp_LRP1.totalb
rs1799986_250kbp_LRP1.totalb
rs2038255_100kbp_FAM177A1.totalb
rs2038255_100kbp_PPP2R3C.totalb
rs2038255_10kbp_FAM177A1.totalb
rs2038255_10kbp_PPP2R3C.totalb
rs2038255_1Mbp_FAM177A1.totalb
rs2038255_1Mbp_PPP2R3C.totalb
rs2038255_250kbp_FAM177A1.totalb
rs2038255_250kbp_PPP2R3C.totalb
rs2038255_3Mbp_FAM177A1.totalb
rs2038255_3Mbp_PPP2R3C.totalb
rs2041733_10kbp_CLEC16A.totalb
rs2143950_100kbp_FAM177A1.totalb
rs2143950_100kbp_PPP2R3C.totalb
rs2143950_10kbp_FAM177A1.totalb
rs2143950_10kbp_PPP2R3C.totalb
rs2143950_1Mbp_FAM177A1.totalb
rs2143950_1Mbp_PPP2R3C.totalb
rs2143950_250kbp_FAM177A1.totalb
rs2143950_250kbp_PPP2R3C.totalb
rs2143950_3Mbp_FAM177A1.totalb
rs2143950_3Mbp_PPP2R3C.totalb
rs2212434_1Mbp_LRRC32.totalb
rs2212434_250kbp_LRRC32.totalb
rs2212434_3Mbp_LRRC32.totalb
rs2218565_100kbp_PRR5L.totalb
rs2218565_10kbp_PRR5L.totalb
rs2218565_1Mbp_PRR5L.totalb
rs2218565_250kbp_PRR5L.totalb
rs2218565_3Mbp_PRR5L.totalb
rs2592555_100kbp_PRR5L.totalb
rs2592555_10kbp_PRR5L.totalb
rs2592555_1Mbp_PRR5L.totalb
rs2592555_250kbp_PRR5L.totalb
rs2592555_3Mbp_PRR5L.totalb
rs3917265_1Mbp_IL1RL1.totalb
rs3917265_3Mbp_IL1RL1.totalb
rs4809219_1Mbp_LIME1.totalb
rs4809219_250kbp_LIME1.totalb
rs4809219_3Mbp_LIME1.totalb
rs6419573_1Mbp_IL1RL1.totalb
rs6419573_250kbp_IL1RL1.totalb
rs6419573_3Mbp_IL1RL1.totalb
rs6602364_100kbp_IL2RA.totalb
rs6602364_1Mbp_IL2RA.totalb
rs6602364_250kbp_IL2RA.totalb
rs6602364_3Mbp_IL2RA.totalb
rs7127307_1Mbp_ETS1.totalb
rs7127307_3Mbp_ETS1.totalb
rs7146581_100kbp_TRAF3.totalb
rs7146581_10kbp_TRAF3.totalb
rs7146581_1Mbp_TRAF3.totalb
rs7146581_250kbp_TRAF3.totalb
rs7146581_3Mbp_TRAF3.totalb

#Generate summary table with significant results.
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
echo -e $snp'\t'$interval'\t'$gene'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>eqtlgen_coloc_results
done < significant_totalb.txt


#Divide the results into one directory for SNPs, each.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_int in 10kbp 100kbp 250kbp 1Mbp 3Mbp  
do
cd $analysis
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
mkdir -p chr${chrom}_${my_snps}/${my_int}
cp -r ${my_snps}*${my_int}* chr${chrom}_${my_snps}/${my_int}
done
done

#Investigation of eQTL pvalues in secondary signals.
for a in rs6419573*3Mbp_*eqtl      
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs3917265*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done


for a in rs6827756*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done


for a in rs13152362*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs12188917*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs4705962*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs2592555*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs2218565*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs12295535*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs61813875*3Mbp_*eqtl 
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

for a in rs7512552*3Mbp_*eqtl
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f4)
echo $gene	$pval >>${snp}_3Mbp_eqtlgen_pvals
done

#Fix spacing to tabs
for b in *_eqtlgen_pvals
do
cat $b | weird >temp
mv temp $b
done

#Find any index SNPs with very low corresponding value in eQTLGen (lower than 10e-50).
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*eqtl
do
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $my_snps | cut -d"_" -f1)
pval=$(grep $snp $my_file | cut -d$'\t' -f4)
echo $snp $gene	$pval >>all_pvals_3Mbp.txt
done
done

sort -k3 -g all_pvals_3Mbp.txt >all_pvals_3Mbp_sorted.txt

#Plot histogram of eqtlgen pvalues for SNPs in the 3 Mbp interval.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*eqtl
do
Rscript $scripts/plot_log_pval.R $my_file ${my_file}_pval.pdf
done
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $my_file | cut -d"_" -f1)
echo $snp $gene	$pval >>${snp}.3Mbp.PPH4
done
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*1Mbp_*totalb
do
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $my_file | cut -d"_" -f1)
echo $snp $gene	$pval >>${snp}.1Mbp.PPH4
done
done

#Fix weird spacing and sort
for b in *.PPH4
do
cat $b | sort -k3 -g | weird >temp
mv temp $b
done

#Prepare PPH4 input for plotting with gassocplot - 1Mbp.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*1Mbp_*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $1}')
snp_chrom=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $2}')
snp_pos=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_1Mbp_eqtlgen.gassocplot2
done
done

#Prepare PPH4 input for plotting with gassocplot - 3Mbp.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $1}')
snp_chrom=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $2}')
snp_pos=$(grep -w $my_snp ${my_file%.totalb}.eqtl | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_3Mbp_eqtlgen.gassocplot2
done
done


#Produce Gassocplots. If for any given interval, we only have Eqtlgen results, we plot a single plot,
#however, if we have both eQTlgen and Sun(2018) then we plot a stack plot.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
if [ -e $sun/${my_snps}_3Mbp_sun.gassocplot ]
then
    Rscript $scripts/plot_stack_gassoc.R ${my_snps}_3Mbp_eqtlgen.gassocplot $sun/${my_snps}_3Mbp_sun.gassocplot \
    ${my_snps}_3Mbp.gwas eQTL pQTL ${my_snps}_3Mbp_stack_gassocplot.pdf
else
    Rscript $scripts/plot_single_gassoc.R ${my_snps}_3Mbp_eqtlgen.gassocplot ${my_snps}_3Mbp.gwas ${my_snps}_3Mbp_single_gassocplot.pdf
fi
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
if [ -e $sun/${my_snps}_1Mbp_sun.gassocplot ]
then
    Rscript $scripts/plot_stack_gassoc.R ${my_snps}_1Mbp_eqtlgen.gassocplot $sun/${my_snps}_1Mbp_sun.gassocplot \
    ${my_snps}_1Mbp.gwas eQTL pQTL ${my_snps}_1Mbp_stack_gassocplot.pdf
else
    Rscript $scripts/plot_single_gassoc.R ${my_snps}_1Mbp_eqtlgen.gassocplot ${my_snps}_1Mbp.gwas ${my_snps}_1Mbp_single_gassocplot.pdf
fi
done

#Generate a list of SNPs closest to the middle of the gene.
python $scripts/middle_of_transcript.py $sun/gencode.v19.annotation.gtf >gencode.v19.annotation.middle
sed -i 's/"//g' gencode.v19.annotation.middle

#Do the same for Ensembl gene names.
python $scripts/middle_of_transcript_ensembl.py $sun/gencode.v19.annotation.gtf >gencode.v19.annotation.ensembl
sed -i 's/"//g' gencode.v19.annotation.ensembl

#Retrieve Ensembl gene names.
Rscript $scripts/get_ensembl_gene_names.R 
#Combine Ensembl and Gencode gene names into one file.
cat gencode.v19.annotation.ensembl.names | awk -v OFS="\t" '{print $3, $2}' | tail -n +2 > gencode_ensembl.v19.annotation.middle
cat gencode.v19.annotation.middle >> gencode_ensembl.v19.annotation.middle
#Remove duplicates.
cat gencode_ensembl.v19.annotation.middle | sort | uniq >temp
mv temp gencode_ensembl.v19.annotation.middle
#Sort the original gassocplot to find the highest pvalue per file (the true hit displayed with the real rsid)

#Prepare input for gassocplot. Modify Gassocplot input files to plot PPH4 values for all the genes in the plot. For the ones that do not meet the significance threshold, plot them in the middle of the gene.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_1Mbp_eqtlgen.gassocplot2 | tail -1 > ${my_snps}_1Mbp_eqtlgen.top
	python $scripts/generate_gassoc_input2.py \
	${my_snps}_1Mbp_eqtlgen.top \
	${my_snps}.1Mbp.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_1Mbp_eqtlgen2
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_3Mbp_eqtlgen.gassocplot2 | tail -1 > ${my_snps}_3Mbp_eqtlgen.top
	python $scripts/generate_gassoc_input2.py \
	${my_snps}_3Mbp_eqtlgen.top \
	${my_snps}.3Mbp.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_3Mbp_eqtlgen2
done


#Plot Sun2018 and eqtlgen in one figure, if possible.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_3Mbp_sun2 ]
then
    Rscript $scripts/plot_stack_gassoc2.R ${my_snps}_3Mbp_eqtlgen2 $sun/${my_snps}_3Mbp_sun2 \
	eQTL pQTL ${my_snps}_3Mbp_stack_gassocplot2.pdf
else
    Rscript $scripts/plot_single_gassoc2.R ${my_snps}_3Mbp_eqtlgen2 ${my_snps}_3Mbp_single_gassocplot2.pdf
fi
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
if [ -e $sun/${my_snps}_1Mbp_sun2 ]
then
    Rscript $scripts/plot_stack_gassoc2.R ${my_snps}_1Mbp_eqtlgen2 $sun/${my_snps}_1Mbp_sun2 \
	eQTL pQTL ${my_snps}_1Mbp_stack_gassocplot2.pdf
else
    Rscript $scripts/plot_single_gassoc2.R ${my_snps}_1Mbp_eqtlgen2 ${my_snps}_1Mbp_single_gassocplot2.pdf
fi
done

#Check if any SNPs in bigld interval is amongst trans-eQTLs from eqtlgen.
#Extract a list of SNPs per each locus to be looked up. 
FINEMAP_ANALYSIS=$HOME/analysis/bayesian_fm/finemap/1k/euro
for locus in $FINEMAP_ANALYSIS/chr*/bigld
do
my_filename=$(basename $locus/*.z | cut -d"." -f1)
cat $locus/*.z | cut -d" " -f1 >$analysis/transeqtl/${my_filename}.bigld
done


cd $analysis/transeqtl
for a in *bigld
do
while read line
do
set -- $line
rsid=$1
echo $rsid
grep -w $rsid $eqtlgen/trans-eQTLs_full_20180905.txt >${a%.bigld}.hits
done < $a
done

#Quicker route to achieve the same.
#All SNPs in BIGLD results
for a in *bigld; do cat $a >>temp; done
cat temp | sort | uniq >all_snps.bigld

#All 10k-odd SNPs in the eqtlgen analysis.
cat $eqtlgen/trans-eQTLs_full_20180905.txt | cut -d" " -f2 | sort | uniq  >$eqtlgen/trans-eQTLs_full_20180905_unique.snps
#Found 133 SNPs featured in the BIgLD interval, as well as amongst the trans-eQTL SNPs in eqtlgen.
#First, find which locus each SNP belongs to.
for a in *.txt *.bigld
do
sort $a >temp
mv temp $a
done

for a in *bigld
do
comm -12 $a bigld_snps_eqtlgen.txt >${a%.bigld}_eqtlgen.bigld
done

for a in rs*_eqtlgen.bigld
do
while read line
do
set -- $line
rsid=$1
echo $rsid
grep -w $rsid $eqtlgen/trans-eQTLs_full_20180905.txt >>${a%.bigld}.hits
done < $a
done

#Sort by pvalues.
for a in *hits
do
sort -k1 -g $a >${a}_sorted
done

for a in *hits_sorted
do
echo $a >>overview_transeqtlgen
head $a >>overview_transeqtlgen
done

#New edition of figures with gassocplot - this time, dot in the middle of gene, for all genes and labelled with the gene name.
#Missing datapoint not plotted.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_1Mbp_eqtlgen.gassocplot2 | tail -1 > ${my_snps}_1Mbp_eqtlgen.top
	python $scripts/generate_gassoc_input3.py \
	${my_snps}_1Mbp_eqtlgen.top \
	${my_snps}.1Mbp.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_1Mbp_eqtlgen3
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_3Mbp_eqtlgen.gassocplot2 | tail -1 > ${my_snps}_3Mbp_eqtlgen.top
	python $scripts/generate_gassoc_input3.py \
	${my_snps}_3Mbp_eqtlgen.top \
	${my_snps}.3Mbp.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_3Mbp_eqtlgen3
done

#Plot Sun2018 and eqtlgen in one figure, if possible.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	echo $my_snps
if [ -e $sun/${my_snps}_3Mbp_sun2 ]
then
    Rscript $scripts/plot_stack_gassoc2.R ${my_snps}_3Mbp_eqtlgen3 $sun/${my_snps}_3Mbp_sun3 \
	eQTL pQTL ${my_snps}_3Mbp_stack_gassocplot3.pdf
else
    Rscript $scripts/plot_single_gassoc2.R ${my_snps}_3Mbp_eqtlgen3 ${my_snps}_3Mbp_single_gassocplot3.pdf
fi
done

for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
if [ -e $sun/${my_snps}_1Mbp_sun2 ]
then
    Rscript $scripts/plot_stack_gassoc2.R ${my_snps}_1Mbp_eqtlgen3 $sun/${my_snps}_1Mbp_sun3 \
	eQTL pQTL ${my_snps}_1Mbp_stack_gassocplot3.pdf
else
    Rscript $scripts/plot_single_gassoc2.R ${my_snps}_1Mbp_eqtlgen3 ${my_snps}_1Mbp_single_gassocplot3.pdf
fi
done
