HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/colocalisation/eqtl_catalogue/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue
onek=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes

cd $analysis

###Remember to use R ver 4.0 to be able to use coloc susie

#For each interval, create a list of SNP rsids present within.
cd GWAS_intervals
for my_gwas in *.gwas
do
rsid=$(echo $my_gwas | cut -d"_" -f1)
echo $rsid
tail -n +2 $my_gwas | cut -f12 > ${rsid}_starting_SNPs.txt
done

cd ..

#We need to create signed LD matrix for all the loci with secondary signal first.
#Starting out with 2 SNPs.

while read line
do
    set -- $line
    plink --r --bfile $onek/ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono \
    --extract GWAS_intervals/${4}_starting_SNPs.txt --matrix --out $4_ld --write-snplist
done < interval_r2_0.2_1k_select.bed 

#Run coloc-susie analysis on loci of choice using Open GWAS datasets.
Rscript --vanilla $scripts/ieugwasr_coloc_susie.R

#Run coloc-susie using eQTL-Catalogue.
for my_r in rs6419573 rs61813875
do
	echo $my_r
	#Rscript --vanilla $scripts/eQTL_Catalogue_coloc_susie_local.R $my_r
	qsub -v my_rsid=$my_r $scripts/sub_run_eqtl_local_susie.sh
done

