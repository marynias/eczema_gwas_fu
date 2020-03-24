#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/paintor
PAINTOR=$HOME/bin/PAINTOR_V3.0
PAINTOR_ANALYSIS=$HOME/analysis/bayesian_fm/paintor/1k/published
analysis=$HOME/analysis/bayesian_fm/RefPanel/1kGenomes
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
#Submit test Paintor run
cd $PAINTOR_ANALYSIS
#qsub $PAINTOR/sub_paintor_test.sh

#Use the locus files and LD files generated in 1kGenomes_ref_panel.sh
#Run Paintor without functional annotations

function paintor_mcmc 
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $PAINTOR_ANALYSIS
mkdir -p chr${chrom}.${snp}.${interval}.mcmc
cd chr${chrom}.${snp}.${interval}.mcmc
cp $analysis/chr${chrom}.${snp}.${interval}.processed ./chr${chrom}.${snp}.${interval}
cp $analysis/chr${chrom}.${snp}.${interval}.ld ./chr${chrom}.${snp}.${interval}.ld
echo "chr${chrom}.${snp}.${interval}" >chr${chrom}.${snp}.${interval}.txt

#Create dummy annotations.
no_snps=$(wc -l chr${chrom}.${snp}.${interval} | cut -d " " -f1)
echo "dummy" >chr${chrom}.${snp}.${interval}.annotations
for (( c=1; c<=${no_snps}; c++)) ; do echo "0" ; done >>chr${chrom}.${snp}.${interval}.annotations
qsub -v input=chr${chrom}.${snp}.${interval}.txt,max_causal=5 $scripts/sub_paintor_mcmc.sh
}


function paintor_exact
{
snp=$1
interval=$2
chrom=$(grep $snp $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $PAINTOR_ANALYSIS
mkdir -p chr${chrom}.${snp}.${interval}.exact
cd chr${chrom}.${snp}.${interval}.exact
cp $analysis/chr${chrom}.${snp}.${interval}.processed ./chr${chrom}.${snp}.${interval}
cp $analysis/chr${chrom}.${snp}.${interval}.ld ./chr${chrom}.${snp}.${interval}.ld
echo "chr${chrom}.${snp}.${interval}" >chr${chrom}.${snp}.${interval}.txt

#Create dummy annotations.
no_snps=$(wc -l chr${chrom}.${snp}.${interval} | cut -d " " -f1)
echo "dummy" >chr${chrom}.${snp}.${interval}.annotations
for (( c=1; c<=${no_snps}; c++)) ; do echo "0" ; done >>chr${chrom}.${snp}.${interval}.annotations
qsub -v input=chr${chrom}.${snp}.${interval}.txt,max_causal=2 $scripts/sub_paintor_exact.sh
}

paintor_mcmc rs112111458 1500000 
paintor_mcmc rs112111458 250000 
paintor_mcmc rs112111458 140000 
paintor_mcmc rs112111458 50000 

paintor_exact rs112111458 1500000
paintor_exact rs112111458 250000
paintor_exact rs112111458 140000
paintor_exact rs112111458 50000

paintor_mcmc rs61813875 1500000
paintor_mcmc rs61813875 250000
paintor_mcmc rs61813875 140000
paintor_mcmc rs61813875 50000

paintor_exact rs61813875 1500000
paintor_exact rs61813875 250000
paintor_exact rs61813875 140000
paintor_exact rs61813875 50000

paintor_mcmc rs2212434 1500000
paintor_mcmc rs2212434 250000
paintor_mcmc rs2212434 140000
paintor_mcmc rs2212434 50000

paintor_exact rs2212434 1500000
paintor_exact rs2212434 250000
paintor_exact rs2212434 140000
paintor_exact rs2212434 50000

#Plot the results with LocusZoom.
#Convert from space-delimited to tab-delimited tables
#Reverse probabilities for the data to be suitable for plotting with GWAS software
for a in chr*mcmc
do
cd $a
cat ${a%.mcmc}.results | tr ' ' '\t' | awk -v OFS="\t" '{(NR==1)?$(NF+1)="pval":$(NF+1)=1-$7; print}' | sort -k7 -g -r >${a%.mcmc}.results.pval
my_snp=$(echo $a | cut -d"." -f2)
my_int=$(echo $a | cut -d"." -f3)
my_int=$(expr $my_int / 1000)
qsub -v input_file=${a%.mcmc}.results.pval,snp_id="rsid",pval="Posterior_Prob",my_ref=$my_snp,my_flank=${my_int}kb,my_prefix=${my_snp}_${my_int}kbp $utils/sub_locus_zoom.sh
cd ../
done

for a in chr*exact
do
cd $a
cat ${a%.exact}.results | tr ' ' '\t' | awk -v OFS="\t" '{(NR==1)?$(NF+1)="pval":$(NF+1)=1-$7; print}' sort -k7 -g -r >${a%.exact}.results.pval
my_snp=$(echo $a | cut -d"." -f2)
my_int=$(echo $a | cut -d"." -f3)
my_int=$(expr $my_int / 1000)
qsub -v input_file=${a%.exact}.results.pval,snp_id="rsid",pval="Posterior_Prob",my_ref=$my_snp,my_flank=${my_int}kb,my_prefix=${my_snp}_${my_int}kbp $utils/sub_locus_zoom.sh
cd ../
done
