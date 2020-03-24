##!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas/twinsuk
scripts=$HOME/bin/eczema_gwas_fu/colocalization/twas
gwas=$HOME/data/gwas/paternoster2015
my_gwas=results.euro.pval.1k
twas=$HOME/bin/fusion_twas-master
utils=$HOME/bin/eczema_gwas_fu/utils
twinsuk=$HOME/analysis/twinsuk
data_manipulation=$HOME/analysis/annotation/data_manipulation

#Prepare TWAS-formated input from TwinsUK - lcl.
#Convert Gen imput to Plink format.
cd $analysis/lcl 
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
plink --gen $twinsuk/lcl/${snp}_3Mbp_all.gen --sample $twinsuk/lcl/${snp}_3Mbp_all.sample --oxford-single-chr $chrom --make-bed --out ${snp}_3Mbp_all
plink --gen $twinsuk/lcl/${snp}_3Mbp_eczema.gen --sample $twinsuk/lcl/${snp}_3Mbp_eczema.sample --oxford-single-chr $chrom  --make-bed --out ${snp}_3Mbp_eczema
plink --gen $twinsuk/lcl/${snp}_3Mbp_noneczema.gen --sample $twinsuk/lcl/${snp}_3Mbp_noneczema.sample --oxford-single-chr $chrom --make-bed --out ${snp}_3Mbp_noneczema
done

#Add gene expression to each gene.
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
python $scripts/generate_twas_input_twinsuk.py $twinsuk/lcl/ensemble_ids_unique_ensembl_names \
$twinsuk/lcl/${snp}_3Mbp_snps_headers_all.rpkm \
${snp}_3Mbp_all.fam \
${snp}_3Mbp_all
python $scripts/generate_twas_input_twinsuk.py $twinsuk/lcl/ensemble_ids_unique_ensembl_names \
$twinsuk/lcl/${snp}_3Mbp_snps_headers_eczema.rpkm \
${snp}_3Mbp_eczema.fam \
${snp}_3Mbp_eczema
python $scripts/generate_twas_input_twinsuk.py $twinsuk/lcl/ensemble_ids_unique_ensembl_names \
$twinsuk/lcl/${snp}_3Mbp_snps_headers_noneczema.rpkm \
${snp}_3Mbp_noneczema.fam \
${snp}_3Mbp_noneczema
done

#Generate matching bed and bim files
for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_all_*.fam
do
ln -s ${file%_3Mbp_all_*.fam}_3Mbp_all.bed ${file%.fam}.bed 
ln -s ${file%_3Mbp_all_*.fam}_3Mbp_all.bim  ${file%.fam}.bim
done
done

for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_eczema_*.fam
do
ln -s ${file%_3Mbp_eczema_*.fam}_3Mbp_eczema.bed ${file%.fam}.bed 
ln -s ${file%_3Mbp_eczema_*.fam}_3Mbp_eczema.bim  ${file%.fam}.bim
done
done

for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_noneczema_*.fam
do
ln -s ${file%_3Mbp_noneczema_*.fam}_3Mbp_noneczema.bed ${file%.fam}.bed 
ln -s ${file%_3Mbp_noneczema_*.fam}_3Mbp_noneczema.bim  ${file%.fam}.bim
done
done

#Important, to make GEMMA output visible to script below:
ln -s ./ output
PATH=$PATH:$twas


for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_all_*.bed
do
    Jobs=$(qstat | grep "Q" | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
       	Jobs=$(qstat | wc -l)
    done
my_in=${file%.bed}
my_temp=${file%.bed}.temp
my_out=${file%.bed}.weights
qsub -v input=$my_in,temp=$my_temp,output=$my_out $scripts/sub_generate_weights.sh
done
done

for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_eczema_*.bed
do
    Jobs=$(qstat | grep "Q" | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
       	Jobs=$(qstat | wc -l)
    done
my_in=${file%.bed}
my_temp=${file%.bed}.temp
my_out=${file%.bed}.weights
qsub -v input=$my_in,temp=$my_temp,output=$my_out $scripts/sub_generate_weights.sh
done
done

for snp in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for file in ${snp}_3Mbp_noneczema_*.bed
do
    Jobs=$(qstat | grep "Q" | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
       	Jobs=$(qstat | wc -l)
    done
my_in=${file%.bed}
my_temp=${file%.bed}.temp
my_out=${file%.bed}.weights
qsub -v input=$my_in,temp=$my_temp,output=$my_out $scripts/sub_generate_weights.sh
done
done

#Create a folder for each subgroup and move files there.
mkdir eczema noneczema all 
mkdir eczema_out noneczema_out all_out 
mv rs*_all_*.RDat ./all/
mv rs*_eczema_*.RDat ./eczema/
mv rs*_noneczema_*.RDat ./noneczema/

#Create POS file.
echo -e PANEL'\t'WGT'\t'ID'\t'CHR'\t'P0'\t'P1'\t'N >$analysis/lcl/eczema.P01.pos
cd $analysis/lcl/eczema
for a in *RDat
do
my_gene=$(ls $a | cut -d"_" -f4 | sed 's/.weights.wgt.RDat//')
#Look up Ensembl ID
my_g=$(grep -w $my_gene $twinsuk/lcl/ensemble_ids_unique_ensembl_names | cut -f2)
my_positions=$(grep -w $my_g $data_manipulation/Homo_sapiens.GRCh37.87.ensembl)
my_chrom=$(echo $my_positions | cut -d" " -f2)
my_start=$(echo $my_positions | cut -d" " -f3)
my_end=$(echo $my_positions | cut -d" " -f4)
echo -e eczema'\t'eczema/$a'\t'$my_gene'\t'$my_chrom'\t'$my_start'\t'$my_end'\t'92 >>$analysis/lcl/eczema.P01.pos
done
(head -1 $analysis/lcl/eczema.P01.pos; tail -n +2 $analysis/lcl/eczema.P01.pos | sort | uniq) >$analysis/lcl/temp
mv $analysis/lcl/temp $analysis/lcl/eczema.P01.pos

echo -e PANEL'\t'WGT'\t'ID'\t'CHR'\t'P0'\t'P1'\t'N >$analysis/lcl/noneczema.P01.pos
cd $analysis/lcl/noneczema
for a in *RDat
do
my_gene=$(ls $a | cut -d"_" -f4 | sed 's/.weights.wgt.RDat//')
#Look up Ensembl ID
my_g=$(grep -w $my_gene $twinsuk/lcl/ensemble_ids_unique_ensembl_names | cut -f2)
my_positions=$(grep -w $my_g $data_manipulation/Homo_sapiens.GRCh37.87.ensembl)
my_chrom=$(echo $my_positions | cut -d" " -f2)
my_start=$(echo $my_positions | cut -d" " -f3)
my_end=$(echo $my_positions | cut -d" " -f4)
echo -e noneczema'\t'noneczema/$a'\t'$my_gene'\t'$my_chrom'\t'$my_start'\t'$my_end'\t'610 >>$analysis/lcl/noneczema.P01.pos
done
(head -1 $analysis/lcl/noneczema.P01.pos; tail -n +2 $analysis/lcl/noneczema.P01.pos | sort | uniq) >$analysis/lcl/temp
mv $analysis/lcl/temp $analysis/lcl/noneczema.P01.pos

echo -e PANEL'\t'WGT'\t'ID'\t'CHR'\t'P0'\t'P1'\t'N >$analysis/lcl/all.P01.pos
cd $analysis/lcl/all
for a in *RDat
do
my_gene=$(ls $a | cut -d"_" -f4 | sed 's/.weights.wgt.RDat//')
#Look up Ensembl ID
my_g=$(grep -w $my_gene $twinsuk/lcl/ensemble_ids_unique_ensembl_names | cut -f2)
my_positions=$(grep -w $my_g $data_manipulation/Homo_sapiens.GRCh37.87.ensembl)
my_chrom=$(echo $my_positions | cut -d" " -f2)
my_start=$(echo $my_positions | cut -d" " -f3)
my_end=$(echo $my_positions | cut -d" " -f4)
echo -e all'\t'all/$a'\t'$my_gene'\t'$my_chrom'\t'$my_start'\t'$my_end'\t'764 >>$analysis/lcl/all.P01.pos
done
(head -1 $analysis/lcl/all.P01.pos; tail -n +2 $analysis/lcl/all.P01.pos | sort | uniq) >$analysis/lcl/temp
mv $analysis/lcl/temp $analysis/lcl/all.P01.pos


cd $analysis/lcl
#Processed GWAS results for chromosome.
my_tissue=( all eczema noneczema )
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v input=$PWD,tissue=$i,chrom=$chr $scripts/sub_run_twas_assoc2.sh
done
done

#Joint/conditional tests and plots
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v input=$PWD,tissue=$i,chrom=$chr $scripts/sub_run_twas_cond2.sh
done
done

#Run coloc
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v input=$PWD,coloc_threshold=0.05,gwasn=103066,tissue=$i,chrom=$chr $scripts/sub_run_twas_coloc2.sh
done
done
########################################################################################
##Results harvesting
#Collect significant results into one file (pre-conditioning)
head -1 $analysis/lcl/all_out/results.euro.pval.twas.chr1.top >twinsuk_lcl.all
for i in "${my_tissue[@]}"
do
for a in $analysis/lcl/${i}_out/*top
do
my_lines=$(wc -l $a | cut -d" " -f1)
if [ $my_lines -gt 1 ]
then
tail -n +2 $a >>twinsuk_lcl.all
fi
done
done
#Collect significant results into one file (after-conditioning)
head -1 $analysis/lcl/all_out/results.euro.pval.twas.chr2.analysis.joint_included.dat >twinsuk_lcl.independent
for i in "${my_tissue[@]}"
do
for a in $analysis/lcl/${i}_out/*included.dat
do
tail -n +2 $a >>twinsuk_lcl.independent
done
done

#Collect significant coloc results
head -1 $analysis/lcl/all_out/results.euro.pval.coloc.chr1 >twinsuk_lcl.coloc.sig0.45
for i in "${my_tissue[@]}"
do
for a in $analysis/lcl/${i}_out/*coloc*
do
cat $a | awk '$25 != "NA" {print $0}' | awk -v OFS="\t" '$25 > 0.45 {print $0}' | tail -n +2 >>twinsuk_lcl.coloc.sig0.45
done
done

