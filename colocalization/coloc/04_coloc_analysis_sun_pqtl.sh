#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/sun_pqtl
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
#Identify genes of interest in the interval around our index SNP and download relevant files.
#Download a GTF file annotation for GChr37.p13 from Gencode.
cd $analysis

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
#Remove chr before chromosome name
sed -i 's/^chr//' gencode.v19.annotation.gtf

#Check if correct gene names used, and if we can find all the gene names from Sun2018 files in there.
cd $HOME/scratch/to_Marie
ls -1 | cut -d"." -f1 >$analysis/sun_gene_list.txt
cd $analysis

while read line
do
set -- $line
g_name=""$1""
if grep --quiet $g_name gencode.v19.annotation.gtf 
then
echo $g_name >>hits
fi
done < sun_gene_list.txt 

#Create BED files with intervals to select for each index SNP.
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-1500000, $3+1500000, $1}' >paternoster_2015_index_snps_sorted_3Mbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-500000, $3+500000, $1}' >paternoster_2015_index_snps_sorted_1Mbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-125000, $3+125000, $1}' >paternoster_2015_index_snps_sorted_250kbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-50000, $3+50000, $1}' >paternoster_2015_index_snps_sorted_100kbp.bed
cat $gwas/paternoster_2015_index_snps_sorted.txt | awk -v OFS="\t" '{print $2, $3-5000, $3+5000, $1}' >paternoster_2015_index_snps_sorted_10kbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_3Mbp.bed
done < paternoster_2015_index_snps_sorted_3Mbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_1Mbp.bed
done < paternoster_2015_index_snps_sorted_1Mbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_250kbp.bed
done < paternoster_2015_index_snps_sorted_250kbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_100kbp.bed
done < paternoster_2015_index_snps_sorted_100kbp.bed

while read line
do
set -- $line
rsid=$4
echo "${line}" >${rsid}_10kbp.bed
done < paternoster_2015_index_snps_sorted_10kbp.bed

#Intersect for gene names
for a in rs*.bed
do
intersectBed -wa -a gencode.v19.annotation.gtf -b $a>>${a%.bed}.gtf
done

#Intersect for SNPs
for a in rs*.bed
do
intersectBed -wa -a $gwas/results.euro.1k.bed -b $a>${a%.bed}.snps.bed
done

#Capitalize alleles in the one file.
for a in ZNF358.10529.19.3.tsv.gz.tab
do
awk -v OFS="\t" '$2 = toupper($2), $3 = toupper($3)' $a >ZNF358.10529.19.3.tsv.gz.tab.cap
done

#Harmonize the effect sizes so that they are relevant to the same allele in both files.
#python $utils/harmonize_beta.py --tab $gwas/results.euro.pval.1k --ref ZNF358.10529.19.3.tsv.gz.tab.cap --header_tab Y --header_ref Y \
#--rsid_tab 12 --rsid_ref 1 --effect_tab 4 --alt_tab 5 --effect_ref 2 --beta_tab 8 --zscore_tab 10 --out results.euro.pval.1k.ea

#Found SNPs with the same ID but multiple effect alleles. For that reason, not proceeding with this type of analysis.

#Subset GWAS results to relevant SNPs in each interval.
for a in rs*snps.bed
do
python $utils/filter_file_by_column.py --tab $a --ref $gwas/results.euro.pval.1k --header_ref Y --rsid_tab 4 --rsid_ref 12 >${a%.snps.bed}.euro.pval
done

#Get the list of all the Gene names to be extracted in a given interval.
for a in rs*.gtf  
do
cat $a | sed 's/^.*\(gene_name .*;\).*$/\1/' | cut -d";" -f1 | sed 's/"//g' | cut -d " " -f2 | sort | uniq >${a%.gtf}.gene_names
done
#Match the alphanumeric part before any delimiters such as "." and "-"
for a in *.gene_names
do
cat $a | cut -d "." -f1 | cut -d "-" -f1 | sort | uniq >${a%}_abbrv
done

#Get a master list of all genes to be analysed for all SNPs (3 Mbp interval)
for a in *3Mbp.gene_names_abbrv
do
cat $a >>all_targets_sun2018_3Mbp
done
cat all_targets_sun2018_3Mbp | sort | uniq >all_targets_sun2018_3Mbp.uniq

while read line
do
set -- $line
g_name=$1
find $HOME/scratch/to_Marie -iname ${g_name}.* -exec cp {} ./ \;
done < all_targets_sun2018_3Mbp.uniq

#Extract relevant SNPs from each file with relevant gene.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_interval in 10kbp 100kbp 250kbp 1Mbp 3Mbp
do
python $scripts/prepare_input_sun2018_pqtl.py ${my_snps}_${my_interval}.gene_names_abbrv ${my_snps}_${my_interval}.euro.pval $my_snps $my_interval 
done
done 

#Delete all empty files (for genes where we don't have pQTL data)
find . -name "*.sun" -size -104c -delete

#Add column with MAF to the file with results.
for a in *euro.pval
do
awk '{print $0, "MAF"}' $a | awk -v OFS="\t" 'NR==1 {print $0}; NR>1 && $7 > 0.5 ? $13 = 1 - $7 : $13 = $7' | sed 's/EFFECT_ALLELE_FREQ$/maf/' | tail -n +2 >temp
mv temp $a
done 

#Run coloc
for my_input in *.sun
do
my_gwas=$(echo $my_input | cut -d"_" -f1,2)
Rscript $scripts/run_coloc_pqtl.R ${my_gwas}.euro.pval $my_input 
done

#Subset to only results with robust support for colocalisation (p > 0.75)
for my_input in *.colocp
do
awk '{if (NR==1) ; else if ($15 > 0.75 && $15 <= 1) print $0}' $my_input >${my_input%.coloc}.sig
done

#Again, remove empty files.
find . -name "*.colocp.sig" -size -1c -delete

#And concatenate significant results into one file, along with overall probablities
for my_input in *sig
do
echo ${my_input%.colocp.sig} >>sun2018.euro.summary.sig
tail -n +2 ${my_input%.colocp.sig}.totalp | tr '\n' '\t' >>sun2018.euro.summary.sig
echo "" >>sun2018.euro.summary.sig
cat $my_input >> sun2018.euro.summary.sig
done

#Grep for nearest genes to the index SNP.
while read line
do
set -- $line
gene=$1
grep -A 1 -B 1 $gene sun2018.euro.summary.sig
done < $gwas/paternoster2015_gene_list

#Generate Locus Plots for selected files and convert to PNG.
qsub -v input_file=$analysis/rs10214237_100kbp_IL7R.colocp,snp_id="snp",pval="SNP.PP.H4",my_ref="rs10214237",my_flank="100kb",my_prefix="rs10214237_100kbp.IL7R.colocp" $utils/sub_locus_zoom.sh 
qsub -v input_file=$analysis/rs17881320_10kbp_STAT3.colocp,snp_id="snp",pval="SNP.PP.H4",my_ref="rs17881320",my_flank="10kb",my_prefix="rs17881320_10kbp.STAT3.colocp" $utils/sub_locus_zoom.sh 

for my_pdf in *IL7R*.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png" 
done 

#Identify the genes with overall significant PP.H4 calculated using p-vals.
for a in *totalp
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done

ls rs4809219_100kbp_TNFRSF6B.totalp
ls rs4809219_1Mbp_TNFRSF6B.totalp
ls rs4809219_250kbp_TNFRSF6B.totalp
ls rs4809219_3Mbp_TNFRSF6B.totalp

#Identify the genes with overall significant PP.H4 calculated using betas and variance of beta.
for a in *totalb
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done

#No difference whatsoever when using beta and its variance
ls rs4809219_100kbp_TNFRSF6B.totalb
ls rs4809219_1Mbp_TNFRSF6B.totalb
ls rs4809219_250kbp_TNFRSF6B.totalb
ls rs4809219_3Mbp_TNFRSF6B.totalb


#Plot with LocusZoom
qsub -v input_file=$analysis/rs4809219_100kbp_TNFRSF6B.colocb,snp_id="snp",pval="SNP.PP.H4",my_ref="rs4809219",my_flank="100kb",my_prefix=r"rs4809219_100kbp_TNFRSF6" $utils/sub_locus_zoom.sh
qsub -v input_file=$analysis/rs4809219_250kbp_TNFRSF6B.colocb,snp_id="snp",pval="SNP.PP.H4",my_ref="rs4809219",my_flank="250kb",my_prefix=r"rs4809219_250kbp_TNFRSF6" $utils/sub_locus_zoom.sh 
qsub -v input_file=$analysis/rs4809219_1Mbp_TNFRSF6B.colocb,snp_id="snp",pval="SNP.PP.H4",my_ref="rs4809219",my_flank="1000kb",my_prefix=r"rs4809219_1Mbp_TNFRSF6" $utils/sub_locus_zoom.sh 
qsub -v input_file=$analysis/rs4809219_3Mbp_TNFRSF6B.colocb,snp_id="snp",pval="SNP.PP.H4",my_ref="rs4809219",my_flank="3000kb",my_prefix=r"rs4809219_3Mbp_TNFRSF6" $utils/sub_locus_zoom.sh  


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

#Change rs145809981 into rs41293864, as the former has been merged into the latter on dbSNP, and this is the ID under which it is present in the Sun dataset.
for a in rs145809981*.euro.pval
do
sed -i 's/rs145809981/rs41293864/g' $a
done

#Run coloc again - perhaps we will get high PP.H4 this time around.
for my_input in rs145809981*.sun
do
my_gwas=$(echo $my_input | cut -d"_" -f1,2)
Rscript $scripts/run_coloc_pqtl.R ${my_gwas}.euro.pval $my_input 
done

#Does not change the results.

##Investigation of pvalues in secondary signals.
for a in rs6419573*3Mbp_*sun      
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs3917265*3Mbp_*sun   
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs6827756*3Mbp_*sun    
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs13152362*3Mbp_*sun   
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs12188917*3Mbp_*sun   
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs4705962*3Mbp_*sun    
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs2592555*3Mbp_*sun  
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs2218565*3Mbp_*sun   
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs12295535*3Mbp_*sun 
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs61813875*3Mbp_*sun 
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

for a in rs7512552*3Mbp_*sun
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $gene $pval $snp >>${snp}_3Mbp_sun_pvals
done

#Fix spacing to tabs
for b in *_sun_pvals
do
cat $b | weird >temp
mv temp $b
done


#Find any index SNPs with very low corresponding value in eQTLGen (lower than 10e-50).
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for a in ${my_snps}*3Mbp_*sun
do
gene=$(echo $a | cut -d"_" -f3 | cut -d"." -f1)
snp=$(echo $a | cut -d"_" -f1)
pval=$(grep $snp $a | cut -d$'\t' -f7)
echo $snp $gene	$pval >>all_pvals_3Mbp_sun.txt
done
done

sort -k3 -g all_pvals_3Mbp_sun.txt | weird >all_pvals_3Mbp_sun_sorted.txt

#Plot histogram of eqtlgen pvalues for SNPs in the 3 Mbp interval.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*sun
do
Rscript $scripts/plot_log_pval_sun.R $my_file ${my_file}_pval.pdf
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


for my_snps in rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2038255 rs6602364  rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs7512552 rs12730935 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*1Mbp_*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.sun | awk -v OFS="\t" '{print $1}')
snp_chrom=$(grep -w $snp_id ${my_snps}_1Mbp.euro.pval | head -1 | awk -v OFS="\t" '{print $2}')
snp_pos=$(grep -w $snp_id ${my_snps}_1Mbp.euro.pval | head -1 | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_1Mbp_sun.gassocplot2
done
done

#Prepare PPH4 input for plotting with gassocplot - 3Mbp.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_file in ${my_snps}*3Mbp_*totalb
do
my_snp=$(awk 'NR == 2 {print $1}' ${my_file%totalb}colocp)
snp_id=$(grep -w $my_snp ${my_file%.totalb}.sun | awk -v OFS="\t" '{print $1}')
snp_chrom=$(grep -w $snp_id ${my_snps}_3Mbp.euro.pval | head -1 | awk -v OFS="\t" '{print $2}')
snp_pos=$(grep -w $snp_id ${my_snps}_3Mbp.euro.pval | head -1 | awk -v OFS="\t" '{print $3}')
pval=$(awk 'NR == 7 {print $2}' $my_file)
gene=$(echo $my_file | cut -d"_" -f3 | cut -d"." -f1)
echo -e $snp_id'\t'$snp_chrom'\t'$snp_pos'\t'$pval'\t'$gene >>${my_snps}_3Mbp_sun.gassocplot2
done
done

#Generate a list of SNPs closest to the middle of the gene.
python $scripts/middle_of_transcript.py $sun/gencode.v19.annotation.gtf >gencode.v19.annotation.middle
sed -i 's/"//g' gencode.v19.annotation.middle
#Do the same but using Ensembl IDs
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
#Sort the original gassocplot to find the highest pvalue per file (the true hit displayed with the real rsid.
#Prepare input for gassocplot. Modify Gassocplot input files to plot PPH4 values for all the genes in the plot. For the ones that do not meet the significance threshold, plot them in the middle of the gene.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_1Mbp_sun.gassocplot2 | tail -1 > ${my_snps}_1Mbp_sun.top
	python $scripts/generate_gassoc_input2.py \
	${my_snps}_1Mbp_sun.top \
	${my_snps}.1Mbp.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_1Mbp_sun2
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_3Mbp_sun.gassocplot2 | tail -1 > ${my_snps}_3Mbp_sun.top
	python $scripts/generate_gassoc_input2.py \
	${my_snps}_3Mbp_sun.top \
	${my_snps}.3Mbp.PPH4 \
	gencode.v19.annotation.middle \
	${my_snps}_3Mbp_sun2
done

#New edition of figures with gassocplot - this time, dot in the middle of gene, for all genes and labelled with the gene name.
#Missing datapoint not plotted.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_1Mbp_sun.gassocplot2 | tail -1 > ${my_snps}_1Mbp_sun.top
	python $scripts/generate_gassoc_input3.py \
	${my_snps}_1Mbp_sun.top \
	${my_snps}.1Mbp.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_1Mbp_sun3
done


for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
	sort -g -k4 ${my_snps}_3Mbp_sun.gassocplot2 | tail -1 > ${my_snps}_3Mbp_sun.top
	python $scripts/generate_gassoc_input3.py \
	${my_snps}_3Mbp_sun.top \
	${my_snps}.3Mbp.PPH4 \
	gencode_ensembl.v19.annotation.middle \
	${my_snps}_3Mbp_sun3
done

