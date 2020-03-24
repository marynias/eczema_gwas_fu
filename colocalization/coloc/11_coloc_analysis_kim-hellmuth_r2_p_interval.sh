#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/working/data/Datasets/eQTL/Kim-Hellmuth2017
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
sun=$HOME/analysis/colocalization/coloc/sun_pqtl
cd $analysis/colocalisation

#Get GWAS results
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d" " -f4)
python $utils/subset_interval_chr.py --tab $gwas/results.euro.pval.1k --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_r2_0.2_1k.gwas
done < $HOME/analysis/annotation/data_manipulation/interval_r2_0.2_1k.bed


done

#Substitute the rsids in both GWAS and colocalisation files to common reference dbsnp RSID.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
python $data_manipulation_scripts/snp_json_substitute.py --tab ${snp}_r2_0.2_1k.gwas \
--db $data_manipulation/rsid_synonyms.txt --head Y  --rsid 12 --delim A >${snp}_r2_0.2_1k.dbsnp.gwas 
done

for sample in lps6h lps90 mdp6h mdp90 ctrl
do
python $data_manipulation_scripts/snp_json_substitute.py --tab ../nominal_${sample}.txt \
--db $data_manipulation/rsid_synonyms.txt --head N --rsid 2 --delim B >nominal_${sample}.dbsnp.txt
done

##Substitute Probe IDs with gene names with ENGS.
#Remove non-data lines
tail -n +21 ../A-MEXP-2210.adf.txt >A-MEXP-2210.adf.data.txt 
for sample in lps6h lps90 mdp6h mdp90 ctrl
do
python $scripts/probe_to_eng_kh.py A-MEXP-2210.adf.data.txt nominal_${sample}.dbsnp.txt $data_manipulation/hugo_synonyms_ids2_filtered.ensembl.upper >nominal_${sample}.dbsnp.eng.txt
done

#Add MAF column
for a in *.dbsnp.gwas 
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 

#Select only lines matching genes in the 3Mbp interval around each SNP
for a in $HOME/analysis/cedar/*.gene_ids
do
my_file=$(basename $a)	
snp=$(echo $my_file | cut -d"_" -f1)
for sample in lps6h lps90 mdp6h mdp90 ctrl
do
while read line
do
set -- $line
gene_id=$1
grep -w "$gene_id" nominal_${sample}.dbsnp.eng.txt >>${snp}_${sample}.dbsnp.eng.txt
done < $a
done
done

#Divide into individual files based on HUGO gene name.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in lps6h lps90 mdp6h mdp90 ctrl
do
python $scripts/extract_kim_hellmuth.py ${snp}_${sample}.dbsnp.eng.txt $data_manipulation/hugo_synonyms_ids2_filtered.hugo ${snp}_${sample}_kh.eqtl
done
done

#Run coloc
function run_coloc
{
sample=$1
sample_no=185
for a in ./*${sample}*_kh.eqtl
do
file=$(basename $a)
snp=$(echo $file | cut -d"_" -f2)
my_output=$(echo $file | sed 's/\.eqtl//')
echo $my_output
Rscript $scripts/coloc_kh.R $sample_no ${my_output}_r2_0.2_1k ${snp}_r2_0.2_1k.dbsnp.gwas $a
done
}

for a in lps6h lps90 mdp6h mdp90 ctrl
do
run_coloc $a
done


for sample in lps6h lps90 mdp6h mdp90 ctrl
do
for a in *${sample}*totalp
do
awk 'NR == 7 && $2 > 0.45 {print FILENAME}' $a >>kh_${sample}_significant_totalp.txt
done
done

function summary_table {
my_tissue=$1
while read line
do
set -- $line
file=$1
echo $file
my_file=$(basename $file)
gene=$(echo $my_file | cut -d"_" -f1)
snp=$(echo $my_file | cut -d"_" -f2)
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_p=$(awk 'NR == 7 {print $2}' $file)
echo -e $gene'\t'$snp'\t'$my_tissue'\t'$pval_p'\t'$target_snp'\t'$target_snp_pp'\t' >>kh_${my_tissue}.coloc_results
done < kh_${my_tissue}_significant_totalp.txt
}

for sample in lps6h lps90 mdp6h mdp90 ctrl
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
tissue=$(echo $my_file | cut -d"_" -f3- | sed 's/_r2_0.2_1k.totalp//')
target_snp=$(awk 'NR == 2 {print $1}' ${file%totalp}colocp)
target_snp_pp=$(awk 'NR == 2 {print $15}' ${file%totalp}colocp)
pval_p=$(awk 'NR == 7 {print $2}' ${file})
echo -e $gene'\t'$snp'\t'$tissue'\t'$pval_p'\t'$target_snp'\t'$target_snp_pp >>${my_tissue}.all.coloc_results
done
}

for a in lps6h lps90 mdp6h mdp90 ctrl
do
summary_table_all_results $a
done
