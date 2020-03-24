#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/finemap
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
onek=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
FINEMAP_ANALYSIS=$HOME/analysis/bayesian_fm/finemap/1k/published
utils=$HOME/bin/eczema_gwas_fu/utils

#Submit test finemap run
cd $FINEMAP_ANALYSIS
#qsub $scripts/sub_finemap_test.sh

#Generate input files for toy analysis fine-mapping the filaggrin locus.
#Important: need to provide MAF (minor allele frequency) rather than EAF (effect allele frequency) like in the Z input file, so need to convert that.

#Copy LD file for European 1K Phase3 along with the list of processed SNPs.
onek=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
#Generate input Z file
input_file=chr1.rs61813875.1500000
python $scripts/generate_Z_file.py --tab $gwas/results.published.tsv \
--proces ${onek}/$input_file.processed --out ${input_file}.z \
--chrom 2 --pos 3 --ident 1 --ref 4 --alt 5 --beta 8 --se 9 --eas 6

#Generate MASTER file
#Using 103066 as sample size (European case and control samples in the discovery phase, including 23 and me. Without 23 and me, that would be 40835)
mkdir $input_file
cd $input_file
mv ../chr1.rs61813875.1500000.z ./ 
echo "z;ld;snp;config;log;n_samples" >master
echo "${input_file}.z;${onek}/${input_file}.ld;${input_file}.snp;${input_file}.config;${input_file}.log;103066" >>master

#shotgun stochastic search
qsub $scripts/sub_finemap.sh

#The results are a bit unexpected - different top loci (except for rs61816766) than in the analysis - perhaps program expects floating point LD numbers rather than scientifc notation?
awk '{for (i=1; i<=NF; i++) printf("%.15f ", $i);} {printf("\n")}' ${onek}/${input_file}.ld > ${onek}/${input_file}_float.ld

mkdir ${input_file}_float
cd ${input_file}_float
cp ../$input_file/chr1.rs61813875.1500000.z ./
echo "z;ld;snp;config;log;n_samples" >master
echo "${input_file}.z;${onek}/${input_file}_float.ld;${input_file}.snp;${input_file}.config;${input_file}.log;103066" >>master
qsub $scripts/sub_finemap.sh
#Obtain exactly the same results.


##Now, analysis of the CD207 locus.
function finemap_default
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $FINEMAP_ANALYSIS
input_file=chr${chrom}.$snp.$interval
mkdir $input_file && cd $input_file
python $scripts/generate_Z_file.py --tab $gwas/results.published.tsv \
--proces ${onek}/$input_file.processed --out ${input_file}.z \
--chrom 2 --pos 3 --ident 1 --ref 4 --alt 5 --beta 8 --se 9 --eas 6
echo "z;ld;snp;config;log;n_samples" >master
echo "${input_file}.z;${onek}/${input_file}.ld;${input_file}.snp;${input_file}.config;${input_file}.log;103066" >>master
qsub $scripts/sub_finemap.sh
}
finemap_default rs112111458 1500000
##Now, analysis of the EMSY/C11orf30 locus.
finemap_default rs2212434 1500000

#Analysis using 500 kbp, 280 kbp and 100 kbp intervals.
finemap_default rs61813875 250000
finemap_default rs112111458 250000
finemap_default rs2212434 250000

finemap_default rs61813875 140000
finemap_default rs112111458 140000
finemap_default rs2212434 140000

finemap_default rs61813875 50000
finemap_default rs112111458 50000
finemap_default rs2212434 50000

#Analysis using only 1 causal SNP for CD207 
function finemap_1snp
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $FINEMAP_ANALYSIS
input_file=chr${chrom}.$snp.$interval
mkdir ${input_file}_1snp && cd ${input_file}_1snp
python $scripts/generate_Z_file.py --tab $gwas/results.published.tsv \
--proces ${onek}/$input_file.processed --out ${input_file}.z \
--chrom 2 --pos 3 --ident 1 --ref 4 --alt 5 --beta 8 --se 9 --eas 6
echo "z;ld;snp;config;log;n_samples" >master
echo "${input_file}.z;${onek}/${input_file}.ld;${input_file}.snp;${input_file}.config;${input_file}.log;103066" >>master
qsub $scripts/sub_finemap_1snp.sh
}
finemap_1snp rs112111458 1500000
finemap_1snp rs112111458 250000
finemap_1snp rs112111458 140000
finemap_1snp rs112111458 50000

#Convert from space-delimited to tab-delimited tables
#Reverse probabilities for the data to be suitable for plotting with GWAS software
#Filter data to keep for plotting
for a in chr*0
do
cd $a
cat ${a}.snp | tr ' ' '\t' >${a}.snp.tab 
cat ${a}.config | tr ' ' '\t' >${a}.config.tab 
awk '($12 > 1) {print $0}' ${a}.snp.tab > ${a}.snp.filtered
awk -v OFS="\t" '{(NR==1)?$(NF+1)="pval":$(NF+1)=1-$11 ; print}' ${a}.snp.filtered >${a}.snp.filtered.pval
my_snp=$(echo $a | cut -d"." -f2)
my_int=$(echo $a | cut -d"." -f3)
my_int=$(expr $my_int / 1000)
qsub -v input_file=${a}.snp.filtered.pval,snp_id="rsid",pval="prob",my_ref=$my_snp,my_flank=${my_int}kb,my_prefix=${my_snp}_${my_int}kbp $utils/sub_locus_zoom.sh
cd ../
done

for a in chr*1snp
do
cd $a
cat ${a%_1snp}.snp | tr ' ' '\t' >${a%_1snp}.snp.tab 
cat ${a%_1snp}.config | tr ' ' '\t' >${a%_1snp}.config.tab 
awk '($12 > 1) {print $0}' ${a%_1snp}.snp.tab > ${a%_1snp}.snp.filtered
awk -v OFS="\t" '{(NR==1)?$(NF+1)="pval":$(NF+1)=1-$11 ; print}' ${a%_1snp}.snp.filtered >${a%_1snp}.snp.filtered.pval
cd ../
done

#Template script to generate LocusPlots and submit them
python $utils/create_batch_scripts.py $scripts/locus_zoom_finemap.sh 

