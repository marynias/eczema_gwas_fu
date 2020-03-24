##!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas/cedar
scripts=$HOME/bin/eczema_gwas_fu/colocalization/twas
gwas=$HOME/data/gwas/paternoster2015
my_gwas=results.euro.pval.1k
twas=$HOME/bin/fusion_twas-master
utils=$HOME/bin/eczema_gwas_fu/utils
cedar=$HOME/analysis/cedar
data_manipulation=$HOME/analysis/annotation/data_manipulation

#Prepare TWAS-formated input.
#Convert Gen imput to Plink format.
cd $analysis

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
plink --gen $cedar/${snp}_3Mbp_${sample}_cedar.gen  --sample $cedar/${snp}_3Mbp_${sample}_cedar.sample --oxford-single-chr $chrom --make-bed --out ${snp}_3Mbp_${sample}_cedar
done
done

#Add gene expression to each gene.
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
python $scripts/generate_twas_input_twinsuk.py $cedar/ensemble_ids_unique_ensembl_names \
$cedar/${snp}_3Mbp_${sample}_GE_Corrected4_Covars_PCs_headers.expr \
${snp}_3Mbp_${sample}_cedar.fam \
${snp}_3Mbp_${sample}_cedar
done
done

#Generate matching bed and bim files
for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
for file in ${snp}_3Mbp_${sample}_cedar_*.fam
do
ln -s ${file%_cedar_*.fam}_cedar.bed ${file%.fam}.bed 
ln -s ${file%_cedar_*.fam}_cedar.bim ${file%.fam}.bim
done
done
done


#Important, to make GEMMA output visible to script below:
ln -s ./ output
PATH=$PATH:$twas

for snp in rs7512552 rs61813875 rs12730935 rs10199605 rs4643526 rs112111458 rs3917265 rs6419573 rs1057258 rs11923593 rs7625909 rs1249910 rs6827756 rs13152362 rs10214237 rs12188917 rs4705962 rs6872156 rs145809981 rs41293864 rs4713555 rs6473227 rs7016497 rs6602364 rs2944542 rs4312054 rs2218565 rs2592555 rs12295535 rs2905493 rs10791824 rs2212434 rs7127307 rs1799986 rs2227483 rs2038255 rs2143950 rs7146581 rs2041733 rs17881320 rs12951971 rs11657987 rs2918307 rs77714197 rs4809219
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do	
for file in ${snp}_3Mbp_${sample}_cedar_*.bed
do
    Jobs=$(qstat | grep "Q" | wc -l)
    while [ $Jobs -gt 100 ]
    do
        sleep 10
        printf "."
       	Jobs=$(qstat | wc -l)
    done
my_in=${file%.bed}
my_temp=${file%.bedqstt}.temp
my_out=${file%.bed}.weights
qsub -v input=$my_in,temp=$my_temp,output=$my_out $scripts/sub_generate_weights.sh
done
done
done

#Create a folder for each subgroup and move files there.
mkdir CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
mkdir CD14_out CD15_out CD19_out CD4_out CD8_out IL_out PLA_out RE_out TR_out
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
for a in rs*_${sample}_*.RDat
do
mv $a $sample
done
done

#Create POS file.
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
echo -e PANEL'\t'WGT'\t'ID'\t'CHR'\t'P0'\t'P1'\t'N >$analysis/${sample}.P01.pos
cd $analysis/$sample
for a in *RDat  
do
my_gene=$(ls $a | cut -d"_" -f5 | sed 's/.weights.wgt.RDat//')
#Look up Ensembl ID
if [[ $my_gene =~ ^ENSG.+$ ]] 
then
my_positions=$(grep -w $my_gene $data_manipulation/Homo_sapiens.GRCh37.87.ensembl)
else 
echo "A"
my_g=$(grep -w $my_gene $HOME/analysis/cedar/ensemble_ids_unique_ensembl_names | cut -f2 )
echo "B"
echo $my_gene
echo $my_g
my_positions=$(grep -w $my_g $data_manipulation/Homo_sapiens.GRCh37.87.ensembl)
fi
echo "C"
my_chrom=$(echo $my_positions | cut -d" " -f2)
my_start=$(echo $my_positions | cut -d" " -f3)
my_end=$(echo $my_positions | cut -d" " -f4)
no_inds=$(wc -l $HOME/analysis/cedar/rs7512552_3Mbp_${sample}_cedar.sample.lst | cut -d" " -f1)
echo -e $sample'\t'$sample/$a'\t'$my_gene'\t'$my_chrom'\t'$my_start'\t'$my_end'\t'$no_inds >>$analysis/${sample}.P01.pos
done
(head -1 $analysis/${sample}.P01.pos; tail -n +2 $analysis/${sample}.P01.pos | sort | uniq) >$analysis/temp
mv $analysis/temp $analysis/${sample}.P01.pos
done

cd $analysis
#Processed GWAS results for chromosome.
for chr in {1..22}
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
qsub -v input=$PWD,tissue=$sample,chrom=$chr $scripts/sub_run_twas_assoc2.sh
done
done

#Joint/conditional tests and plots
for chr in {1..22}
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
qsub -v input=$PWD,tissue=$sample,chrom=$chr $scripts/sub_run_twas_cond2.sh
done
done

#Run coloc
for chr in {1..22}
do
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
qsub -v input=$PWD,coloc_threshold=0.05,gwasn=103066,tissue=$sample,chrom=$chr $scripts/sub_run_twas_coloc2.sh
done
done

########################################################################################
##Results harvesting
#Collect significant results into one file (pre-conditioning).
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
head -1 $analysis/CD14_out/results.euro.pval.twas.chr20.top >cedar_${sample}.all
for a in $analysis/${sample}_out/*top
do
my_lines=$(wc -l $a | cut -d" " -f1)
if [ $my_lines -gt 1 ]
then
tail -n +2 $a >>cedar_${sample}.all
fi
done
done

#Collect significant results into one file (after-conditioning)
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
head -1 $analysis/CD4_out/results.euro.pval.twas.chr20.analysis.joint_included.dat >cedar_${sample}.independent
for a in $analysis/${sample}_out/*included.dat
do
tail -n +2 $a >>cedar_${sample}.independent
done
done

#Collect significant coloc results
for sample in CD14 CD15 CD19 CD4 CD8 IL PLA RE TR
do
head -1 $analysis/CD4_out/results.euro.pval.coloc.chr20 >cedar_${sample}.coloc.sig0.45
for a in $analysis/${sample}_out/*coloc*
do
cat $a | awk '$25 != "NA" {print $0}' | awk -v OFS="\t" '$25 > 0.45 {print $0}' | tail -n +2 >>cedar_${sample}.coloc.sig0.45
done
done

