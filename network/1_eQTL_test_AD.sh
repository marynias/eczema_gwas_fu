#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/network/AD
scripts=$HOME/bin/eczema_gwas_fu/network
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
utils=$HOME/bin/eczema_gwas_fu/utils
coloc_scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
cd $analysis

##Using the six tissues in the GTEx dataset.
#Take known AD loci and artificially lower GWAS p-values to 1e-4 5e-8 level and see if we detect any colocalisation with coloc.

#Take known AD loci and artificially lower GWAS AND eQTL p-values to 1e-4 5e-8 level and see if we detect any colocalisation with coloc.

#Take loci with p-values of between 1e-4 and 5e-8 to see if we detect any colocalisation with coloc.
cat $gwas/results.euro.pval.gtex | (head -1; awk -v OFS="\t" '($11 < 1e-4 && $11 > 5e-8) {print}') | sort -k2,2n -k3,3n>results.euro.highpval
cat $gwas/results.euro.pval.gtex | (head -1; awk -v OFS="\t" '($11 < 1e-6 && $11 > 5e-8) {print}') | sort -k2,2n -k3,3n>results.euro.lowpval
cat $gwas/results.euro.pval.gtex | (head -1; awk -v OFS="\t" '($11 < 1e-5 && $11 > 5e-8) {print}') | sort -k2,2n -k3,3n>results.euro.middlepval
cat $gwas/results.euro.pval.gtex | (head -1; awk -v OFS="\t" '($11 < 0.05 && $11 > 1e-4) {print}') | sort -k2,2n -k3,3n>results.euro.lowestpval

tail -n +2 results.euro.lowestpval | awk -v OFS="\t" '{print $2, $3, $3+1, $4, $5, $8, $9, $11, $12}' >results.euro.lowestpval.bed

#Downloaded USCS cytoband annotation in bed format from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
#Annotate with cytobands
sed 's/chr//g' ucsc_hg19_cytoBand >ucsc_hg19_cytoBand.bed
intersectBed -loj -a results.euro.highpval.bed -b ucsc_hg19_cytoBand.bed >results.euro.highpval.cytoband
intersectBed -loj -a results.euro.lowpval.bed -b ucsc_hg19_cytoBand.bed >results.euro.lowpval.cytoband
intersectBed -loj -a results.euro.middlepval.bed -b ucsc_hg19_cytoBand.bed >results.euro.middlepval.cytoband
intersectBed -loj -a results.euro.lowestpval.bed -b ucsc_hg19_cytoBand.bed >results.euro.lowestpval.cytoband
#New field giving precise cytoband location.
cat results.euro.highpval.cytoband | awk -v OFS="\t" '{print $0, $14$17}' >results.euro.highpval.cytoband2
cat results.euro.lowpval.cytoband | awk -v OFS="\t" '{print $0, $14$17}' >results.euro.lowpval.cytoband2
cat results.euro.middlepval.cytoband | awk -v OFS="\t" '{print $0, $14$17}' >results.euro.middlepval.cytoband2
cat results.euro.lowestpval.cytoband | awk -v OFS="\t" '{print $0, $10$13}' > s
#Using middle threshold results in adding in around 28 new loci, while high threshold around 200 new loci. Sticking with the high threshold for now. 
#Using conservative clumping - selecting one index SNP per broadly defined cytoband, e.g. 21q22. Results in 128 loci for further investigation with colocalisation (250 kbp around the index SNP).

#Subset interval.
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d$' ' -f4)
echo $rsid
python $utils/subset_interval_chr.py --tab $gwas/results.euro.pval.gtex --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_250kbp.gwas
done < new_ad_candidates2.bed

while read line
do
echo $line >temp
rsid=$(cat temp | cut -d$' ' -f4)
echo $rsid
python $utils/subset_interval_chr.py --tab $gwas/results.euro.pval.gtex --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_250kbp.gwas
done < lowest_candidates2.bed

#Add MAF column
for a in *gwas
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 


while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_3Mbp.bed
done < new_ad_candidates3.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_3Mbp.bed
done < lowest_candidates3.bed



#Intersect for gene names
for a in rs*.bed 
do
intersectBed -wa -a $HOME/analysis/annotation/data_manipulation/gtf/gencode.v19.annotation.gtf -b $a>>${a%.bed}.gtf
done

#Get the list of all the Gene names to be extracted in a given interval.
for a in rs*.gtf  
do
cat $a | sed 's/^.*\(gene_id .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | cut -d "." -f1 |sort | uniq >${a%.gtf}.gene_names
done

#Again, for the combined file
for a in rs*.gtf
do
my_file=$(basename $a)
my_rsid=$(echo $my_file | cut -d"_" -f1)
cat $a | sed 's/^.*\(gene_id .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | cut -d"." -f1 | sort | uniq | awk -v var="$my_rsid" '{print $0, var}' | weird >>ensembl_250kbp_interval.gene_ids
done

for a in *.gene_names
do
cat $a >>all_targets_250kbp
done

cat all_targets_250kbp | sort | uniq >all_targets_250kbp.uniq


for my_gtex in Cells_EBV-transformed_lymphocytes.allpairs.txt.gz Cells_Transformed_fibroblasts.allpairs.txt.gz \
Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz \
Spleen.allpairs.txt.gz Whole_Blood.allpairs.txt.gz
do
my_t=$(echo $my_gtex | cut -d "." -f1)
qsub -v gtex_data=$my_gtex,ensembl_selection=ensembl_250kbp_interval.gene_ids,gene_names=$sun/gencode.v19.annotation.ensembl.names,gtex_ref=GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt,tissue=$my_t $coloc_scripts/sub_extract_gtex.sh 
done

for my_gtex in Cells_EBV-transformed_lymphocytes.allpairs.txt.gz Cells_Transformed_fibroblasts.allpairs.txt.gz \
Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz \
Spleen.allpairs.txt.gz Whole_Blood.allpairs.txt.gz
do
my_t=$(echo $my_gtex | cut -d "." -f1)
qsub -v gtex_data=$my_gtex,ensembl_selection=ensembl_3Mbp_interval.gene_ids,gene_names=$sun/gencode.v19.annotation.ensembl.names,gtex_ref=GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt,tissue=$my_t $coloc_scripts/sub_extract_gtex.sh 
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
#Rscript $coloc_scripts/coloc_gtex.R $sample_no ${my_output}_250kbp ../${snp}_250kbp.gwas $a
Rscript $coloc_scripts/coloc_gtex.R $sample_no ${my_output}_3Mbp ../${snp}_250kbp.gwas $a
done
}

run_coloc $analysis/Skin_Not_Sun_Exposed_Suprapubic 335
run_coloc $analysis/Skin_Sun_Exposed_Lower_leg 414
run_coloc $analysis/Whole_Blood 369
run_coloc $analysis/Cells_EBV-transformed_lymphocytes 117
run_coloc $analysis/Cells_Transformed_fibroblasts 300
run_coloc $analysis/Spleen 146

function pvals_p {
my_tissue=$1
#for a in $my_tissue/*${my_tissue}_250kbp.totalp
#do
#awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}_250kbp.sigp
done
for a in $my_tissue/*${my_tissue}_3Mbp.totalp
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.75 {print FILENAME}' $a >>${my_tissue}_3Mbp.sigp
done
}

function pvals_b {
my_tissue=$1
#for a in $my_tissue/*${my_tissue}_250kbp.totalb
#do
#awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}_250kbp.sigb
#done
for a in $my_tissue/*${my_tissue}_3Mbp.totalb
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.75 {print FILENAME}' $a >>${my_tissue}_3Mbp.sigb
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

function summary_table {
my_tissue=$1
while read line
do
set -- $line
file=$1
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f2)
snp=$(echo $my_file | cut -d"_" -f1)
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalb}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalb}colocp)
pval_b=$(awk 'NR == 7 {print $2}' $file)
pval_p=$(awk 'NR == 7 {print $2}' ${file%b}p)
echo -e $gene'\t'$snp'\t'$my_tissue'\t'$pval_p'\t'$pval_b'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}_3Mbp.coloc_results
done < ${my_tissue}_3Mbp.sigb
}

for a in Skin_Not_Sun_Exposed_Suprapubic \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes \
Spleen
do
summary_table $a
done

for my_tissue in Skin_Not_Sun_Exposed_Suprapubic \
Cells_Transformed_fibroblasts \
Skin_Sun_Exposed_Lower_leg \
Whole_Blood \
Cells_EBV-transformed_lymphocytes \
Spleen
do
for a in $my_tissue/*${my_tissue}_3Mbp.totalp
do
filename=$(basename $a)
#echo $filename
snp=$(echo $filename | cut -d"_" -f1)
PPH4=$(awk 'NR == 7 {print $2}' $a)
echo -e $snp'\t'$PPH4'\t'$my_tissue >>${my_tissue}_3Mbp.allp
done
done

for h4 in *allp
do
cat $h4 >>pph4_3mbp_all_combined.txt
done


