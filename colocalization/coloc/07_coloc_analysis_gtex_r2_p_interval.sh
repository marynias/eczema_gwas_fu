#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/gtex
ref_panel=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
twinsuk_scripts=$HOME/bin/eczema_gwas_fu/analysis/twinsuk
qtl_data=$HOME/data/eqtl/gtex
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
utils=$HOME/bin/eczema_gwas_fu/utils
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
cd $analysis

#Download reference file with variants from https://gtexportal.org/home/datasets
wget https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz

#Generate a GTEx version of the GWAS results. 
for a in GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt
do
tail -n +2 $a | awk -v OFS="\t" '{print $1, $7, '0', $2, $4, $5}' >gtex.bim
done

python $utils/update_rsid.py --bim gtex.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$gwas/results.euro.pval.gtex

#Subset interval.
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d" " -f4)
python $utils/subset_interval_chr.py --tab $gwas/results.euro.pval.gtex --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_r2_0.2_1k.gwas
done < $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed

#Add MAF column
for a in *gwas
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 
#Prepare eQTL input for coloc
#N samples in different tissues:

#Whole_Blood	369
#Spleen	146 
#Sun_Exposed_Lower_leg	414
#Skin_Not_Sun_Exposed_Suprapubic	335
#Cells_Transformed_fibroblasts	300
#Cells_EBV-transformed_lymphocytes	117

#Filter results for genes of interest in each dataset.
#First, extract gene Ensembl IDS of interest.
for a in $sun/rs*3Mbp*.gtf
do
my_file=$(basename $a)
my_rsid=$(echo $my_file | cut -d"_" -f1)
cat $a | sed 's/^.*\(gene_id .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | cut -d"." -f1 | sort | uniq | awk -v var="$my_rsid" '{print $0, var}' | weird >>ensembl_3Mbp_interval.gene_ids
done

#Temporarily copy the input eQTL files onto cluster.
for my_gtex in Cells_EBV-transformed_lymphocytes.allpairs.txt.gz Cells_Transformed_fibroblasts.allpairs.txt.gz \
Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz \
Spleen.allpairs.txt.gz Whole_Blood.allpairs.txt.gz
do
my_t=$(echo $my_gtex | cut -d "." -f1)
qsub -v gtex_data=$my_gtex,ensembl_selection=ensembl_3Mbp_interval.gene_ids,gene_names=$sun/gencode.v19.annotation.ensembl.names,gtex_ref=GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt,tissue=$my_t $scripts/sub_extract_gtex.sh 
done

#Divide into individual folders per tissue
for a in Skin_Not_Sun_Exposed_Suprapubic \
Spleen \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes
do
mv *${a}* $a
done

#Run coloc
function run_coloc
{
my_path=$1
sample_no=$2
cd $my_path
for a in ./*.gtex
do
file=$(basename $a)
snp=$(echo $file | cut -d"_" -f1)
my_output=$(echo $file | cut -d"_" -f1,2,4- | sed 's/\.gtex//')
Rscript $scripts/coloc_gtex.R $sample_no ${my_output}_r2_0.2_1k ../${snp}_r2_0.2_1k.gwas $a
done
}

run_coloc $analysis/Skin_Not_Sun_Exposed_Suprapubic 335
run_coloc $analysis/Skin_Sun_Exposed_Lower_leg 414
run_coloc $analysis/Whole_Blood 369
run_coloc $analysis/Cells_EBV-transformed_lymphocytes 117
run_coloc $analysis/Cells_Transformed_fibroblasts 300
run_coloc $analysis/Spleen 146

#Check the colocalisation results.
function pvals_p {
my_tissue=$1
for a in $my_tissue/*${my_tissue}_250kbp.totalp
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}.sigp
done
}

function pvals_b {
my_tissue=$1
for a in $my_tissue/*${my_tissue}_250kbp.totalb
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}.sigb
done
}

for a in Skin_Not_Sun_Exposed_Suprapubic \
Spleen \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes
do
pvals_p $a
pvals_b $a
done

#Generate summary table with significant results. Using p-value based results for lymphocytes and spleen, as low number of samples, which often results in NA - not able to calculate colocalisation using beta and its variance.
function summary_table {
my_tissue=$1
while read line
do
set -- $line
file=$1
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
tissue=$(echo $my_file | cut -d"_" -f3- | sed 's/_r2_0.2_1k.totalb//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $file)
pval_p=$(awk 'NR == 7 {print $2}' ${file%b}p)
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.coloc_results
done < ${my_tissue}.sigb
}


for a in Skin_Not_Sun_Exposed_Suprapubic \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood 
do
summary_table $a
done

function summary_table2 {
my_tissue=$1
while read line
do
set -- $line
file=$1
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
tissue=$(echo $my_file | cut -d"_" -f3- | sed 's/_r2_0.2_1k.totalp//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_b=$(awk 'NR == 7 {print $2}' ${file%p}b)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.coloc_results
done < ${my_tissue}.sigp
}

for a in Cells_EBV-transformed_lymphocytes \
Spleen 
do
summary_table2 $a
done

function summary_table_all_results {
my_tissue=$1
echo -e "Gene\tLead_SNP_interval\tTissue\tTOTAL.PP.H4_pval\tTOTAL.PP.H4_beta\tTop_SNP\tTop_SNP.PP.H4" >>${my_tissue}.all.coloc_results
for file in ${my_tissue}/*.totalp
do
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
tissue=$(echo $my_file | cut -d"_" -f3- | sed 's/_r2_0.2_1k.totalp//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_b=$(awk 'NR == 7 {print $2}' ${file%p}b)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.all.coloc_results
done
}

for a in Skin_Not_Sun_Exposed_Suprapubic \
Spleen \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes
do
summary_table_all_results $a
done



#Prepare PPH4 summary file. First beta-based, then p-value based
function prepare_pph4
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in $tissue/${my_snps}_*${tissue}_r2_0.2_1k.totalb
do
my_file=$(basename $file)
pval=$(awk 'NR == 7 {print $2}' $file)
gene=$(echo $my_file | cut -d"_" -f2 | sed 's/_${tissue}_r2_0.2_1k.totalb//')
echo $my_snps $gene	 $pval >>${my_snps}_gtex_${tissue}_r2_0.2_1k.PPH4
done
done
}

for a in Skin_Not_Sun_Exposed_Suprapubic \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood 
do
prepare_pph4 $a
done

function prepare_pph4_2
{
tissue=$1
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in $tissue/${my_snps}_*${tissue}_r2_0.2_1k.totalp
do
my_file=$(basename $file)
pval=$(awk 'NR == 7 {print $2}' $file)
gene=$(echo $my_file | cut -d"_" -f2 | sed 's/_${tissue}_r2_0.2_1k.totalp//')
echo $my_snps $gene	 $pval >>${my_snps}_gtex_${tissue}_r2_0.2_1k.PPH4
done
done
}

for a in Cells_EBV-transformed_lymphocytes \
Spleen 
do
prepare_pph4_2 $a
done

for my_tissue in Skin_Not_Sun_Exposed_Suprapubic \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes \
Spleen
do
for a in $my_tissue/*${my_tissue}_r2_0.2_1k.totalp
do
filename=$(basename $a)
#echo $filename
snp=$(echo $filename | cut -d"_" -f1)
PPH4=$(awk 'NR == 7 {print $2}' $a)
echo -e $snp'\t'$PPH4'\t'$my_tissue >>${my_tissue}_r2_0.2_1k.allp
done
done

for h4 in *allp
do
cat $h4 >>pph4_r2_0.2_1k_all_combined_sig.txt
done