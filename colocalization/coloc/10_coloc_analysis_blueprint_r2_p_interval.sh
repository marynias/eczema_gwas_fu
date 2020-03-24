#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/working/data/Datasets/eQTL/Chen2016
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
cd $analysis/colocalisation

#Prepare gene-based input for Blueprint: mono, neut and tcel.
#First create a BIM reference file for each file.
for a in tcel mono neut
do
cat $analysis/${a}_gene_nor_combat_peer_10_all_summary.txt | cut -d" " -f1-2 | sed 's/[:_]/\t/g' | awk -v OFS='\t' '{print $1, $5, 0, $2, $3, $4 }' > $analysis/colocalisation/$a/${a}_gene_nor_combat_peer_10_all_summary.bim
done

for a in tcel mono neut
do
cat $analysis/colocalisation/$a/${a}_gene_nor_combat_peer_10_all_summary.bim | sort | uniq >temp
mv temp $analysis/colocalisation/$a/${a}_gene_nor_combat_peer_10_all_summary.bim 
done

#Generate a Blueprint version of the GWAS results. 
for a in tcel mono neut
do
python $utils/update_rsid.py --bim $analysis/colocalisation/$a/${a}_gene_nor_combat_peer_10_all_summary.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$analysis/colocalisation/$a/results.euro.pval.blueprint
done

for a in tcel mono neut
do
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d" " -f4)
python $utils/subset_interval_chr.py --tab $analysis/colocalisation/$a/results.euro.pval.blueprint --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output $analysis/colocalisation/$a/${rsid}_r2_0.2_1k.gwas
done < $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed
done

#Add MAF column
for b in tcel mono neut
do
for a in $analysis/colocalisation/$b/*gwas
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 
done

#Prepare eQTL input for coloc
#N samples in different tissues:

#197

#Filter results for genes of interest in each dataset.
#First, extract gene Ensembl IDS of interest.
for a in $sun/rs*3Mbp*.gtf
do
my_file=$(basename $a)
my_rsid=$(echo $my_file | cut -d"_" -f1)
cat $a | sed 's/^.*\(gene_id .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | cut -d"." -f1 | sort | uniq | awk -v var="$my_rsid" '{print $0, var}' | weird >>ensembl_3Mbp_interval.gene_ids
done

cd $analysis/colocalisation
for a in tcel mono neut
do
python $scripts/extract_blueprint.py $analysis/${a}_gene_nor_combat_peer_10_all_summary.txt \
$analysis/colocalisation/ensembl_3Mbp_interval.gene_ids \
$sun/gencode.v19.annotation.ensembl.names $a 
done


#Divide into individual folders per tissue
for a in tcel mono neut
do
mv *${a}.blueprint $a
done

#Run coloc
function run_coloc
{
my_path=$1
sample_no=$2
cd $my_path
for a in ./rs*.blueprint
do
file=$(basename $a)
snp=$(echo $file | cut -d"_" -f1)
my_output=$(echo $file | cut -d"_" -f1,2,4- | sed 's/\.blueprint//')
echo $my_output
Rscript $scripts/coloc_blueprint.R $sample_no ${my_output}_r2_0.2_1k ${snp}_r2_0.2_1k.gwas $a
done
}

for a in tcel mono neut
do
run_coloc $analysis/colocalisation/$a 197
done

function pvals_p {
my_tissue=$1
for a in $my_tissue/*${my_tissue}_r2_0.2_1k.totalp
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}.sigp
done
}

function pvals_b {
my_tissue=$1
for a in $my_tissue/*${my_tissue}_r2_0.2_1k.totalb
do
awk 'NR == 7 && $2 !~ "NA" && $2 > 0.45 {print FILENAME}' $a >>${my_tissue}.sigb
done
}

for a in tcel mono neut
do
pvals_p $a
pvals_b $a
done

#Generate summary table with significant results using beta-based results.
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

for a in tcel mono neut
do
summary_table_all_results $a
done



for a in tcel mono neut
do
summary_table $a
done

#Prepare PPH4 summary file. 
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
echo $my_snps $gene	 $pval >>${my_snps}_blueprint_${tissue}_r2_0.2_1k.PPH4
done
done
}

for a in tcel mono neut
do
prepare_pph4 $a
done
