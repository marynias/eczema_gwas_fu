#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/regfm
scripts=$HOME/bin/eczema_gwas_fu/annotation/regfm
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation

cd $analysis

#Run test regfm job
qsub -v input="example" $scripts/sub_regfm.sh 

#Link PLINK files with 1k reference to correct folder and filenames for regfm to run.
cd /panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
for chrom in {1..22}
do
ln -s $PWD/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim \
$HOME/bin/regfm/BigData/1000Genome/European/ALL_chr${chrom}.bim
ln -s $PWD/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bed \
$HOME/bin/regfm/BigData/1000Genome/European/ALL_chr${chrom}.bed
ln -s $PWD/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.fam \
$HOME/bin/regfm/BigData/1000Genome/European/ALL_chr${chrom}.fam
done

#Generate input files for regfm finemapping. Use SNPs in predefined intervals.
#Define lead SNPs. Exclude duplicate SNP: rs145809981
echo -e "SNP"'\t'"Chr"'\t'"Position"'\t'"P" >$HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt
while read line
do
set -- $line
rsid=$1
chrom=$2
pos=$3
pval=$(grep -w $rsid $gwas/results.euro.pval.1k | cut -f11)
echo -e $rsid'\t'$chrom'\t'$pos'\t'$pval >>$HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt
done < $gwas/paternoster_2015_index_snps_sorted.txt
sed -i '/^rs145809981/d' $HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt

#Use intervals - first 1k_r2_0.2, then 1Mbp and 3Mbp.
a=$data_manipulation/interval_r2_0.2_1k.gwas
cat $a | awk -v OFS="\t" '{print $2, $3-1, $3, $12, $11}' | sed 's/^/chr/g' | tail -n +2 >$HOME/bin/regfm/BedData/paternoster2015_1k/allSNPs.bed
qsub -v input="paternoster2015_1k" $scripts/sub_regfm.sh 

mkdir -p $HOME/bin/regfm/BedData/paternoster2015_1Mbp
mkdir -p $HOME/bin/regfm/BedData/paternoster2015_3Mbp
cp $HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt $HOME/bin/regfm/BedData/paternoster2015_1Mbp
cp $HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt $HOME/bin/regfm/BedData/paternoster2015_3Mbp
a=$data_manipulation/paternoster_2015_index_snps_sorted_1Mbp.gwas
cat $a | awk -v OFS="\t" '{print $2, $3-1, $3, $12, $11}' | sed 's/^/chr/g' | tail -n +2 >$HOME/bin/regfm/BedData/paternoster2015_1Mbp/allSNPs.bed
a=$data_manipulation/paternoster_2015_index_snps_sorted_3Mbp.gwas
cat $a | awk -v OFS="\t" '{print $2, $3-1, $3, $12, $11}' | sed 's/^/chr/g' | tail -n +2 >$HOME/bin/regfm/BedData/paternoster2015_3Mbp/allSNPs.bed
qsub -v input="paternoster2015_1Mbp" $scripts/sub_regfm.sh 
qsub -v input="paternoster2015_3Mbp" $scripts/sub_regfm.sh 

#Does not work. For testing, use only one locus and the whole GWAS dataset with results.
head -2 $HOME/bin/regfm/BedData/paternoster2015_1k/LeadSNPs.txt >$HOME/bin/regfm/BedData/paternoster2015_1locus_test/LeadSNPs.txt
awk -v OFS="\t" '{print $2, $3-1, $3, $12, $11}' $gwas/results.euro.pval.1k | sed 's/^/chr/g' | tail -n +2 >$HOME/bin/regfm/BedData/paternoster2015_1locus_test/allSNPs.bed
qsub -v input="paternoster2015_1locus_test" $scripts/sub_regfm.sh 