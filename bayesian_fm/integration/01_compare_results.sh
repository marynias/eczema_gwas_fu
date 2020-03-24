#!/bin/bash
scripts=/panfs/panasas01/sscm/qh18484/bin/eczema_gwas_fu/bayesian_fm/integration
paintor=/panfs/panasas01/sscm/qh18484/bin/PAINTOR_V3.0/PAINTOR_Utilities
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
temp=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/integration

cd $temp

#Generate position map for all SNPs analysed - both in the UKBiobank and 1k. Remove duplicates.
cat $gwas/results.euro.pval.1k | awk -v OFS="\t" '{print $12, $2, $3}'  >onek.map
cat $gwas/results.euro.pval.ukbiobank | awk -v OFS="\t" '{print $12, $2, $3}' >ukbiobank.map
cat onek.map | tail -n +2 >combined.map
cat ukbiobank.map | tail -n +2  >>combined.map
sort combined.map | uniq >combined_uniq.map

#Compare the results of JAM, Paintor and Finemap.
#First, pick up high PP targets to compare results across software runs.
finemap=$HOME/analysis/bayesian_fm/finemap/1k/euro/copy
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $finemap/chr${chrom}.${my_snps}/${my_int}
#The last two arguments are LogBF threshold and posterior probability threshold.
python $scripts/grab_top_hits_finemap.py $chrom $my_snps $my_int 2 0.1 >>$temp/chr$chrom.${my_snps}.rsid
cd $temp
done
done

#For now adding loci from extended credible set is disabled.
jam=$HOME/analysis/bayesian_fm/jam/ukbiobank/30k
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
#The last two arguments are BF threshold and posterior probability threshold.
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $jam/chr${chrom}.${my_snps}/${my_int}
python $scripts/grab_top_hits_jam.py $chrom $my_snps $my_int 100 0.1 >>$temp/chr$chrom.${my_snps}.rsid
cd $temp
done
done

paintor=$HOME/analysis/bayesian_fm/paintor/1k/euro
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
#The last two arguments are BF threshold and posterior probability threshold.
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $paintor/chr${chrom}.${my_snps}/${my_int}.exact
python $scripts/grab_top_hits_paintor.py $chrom $my_snps $my_int 100 0.1 >>$temp/chr$chrom.${my_snps}.rsid
cd $paintor/chr${chrom}.${my_snps}/${my_int}.mcmc
python $scripts/grab_top_hits_paintor.py $chrom $my_snps $my_int 100 0.1 >>$temp/chr$chrom.${my_snps}.rsid
cd $temp
done
done

#Remove duplicate hits.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cat $temp/chr$chrom.${my_snps}.rsid | sort | uniq > $temp/chr$chrom.${my_snps}_unique.rsid
done


#Generate a matrix with PP values from selected SNPs for clustering. First generate a dir for the results.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
mkdir -p $temp/chr$chrom.${my_snps}_clustering
done

#Generate a matrix with PP values from selected SNPs for clustering.
#Finemap
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $finemap/chr${chrom}.${my_snps}/${my_int}
python $scripts/generate_input_finemap.py $chrom $my_snps $my_int \
$temp/chr$chrom.${my_snps}_unique.rsid \
$temp/combined_uniq.map \
$temp/chr$chrom.${my_snps}_clustering/finemap_chr$chrom.${my_snps}.${my_int}
cd $temp
done
done


#JAM
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $jam/chr${chrom}.${my_snps}/${my_int}
python $scripts/generate_input_jam.py $chrom $my_snps $my_int \
$temp/chr$chrom.${my_snps}_unique.rsid \
$temp/combined_uniq.map \
$temp/chr$chrom.${my_snps}_clustering/jam_chr$chrom.${my_snps}.${my_int}
cd $temp
done
done

#Paintor
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
for my_int in 5000 50000 250000 500000 1500000 gpart bigld
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $paintor/chr${chrom}.${my_snps}/${my_int}.exact
python $scripts/generate_input_paintor.py $chrom $my_snps $my_int \
$temp/chr$chrom.${my_snps}_unique.rsid \
$temp/combined_uniq.map \
$temp/chr$chrom.${my_snps}_clustering/paintor_chr$chrom.${my_snps}.${my_int}.exact
cd $paintor/chr${chrom}.${my_snps}/${my_int}.mcmc
python $scripts/generate_input_paintor.py $chrom $my_snps $my_int \
$temp/chr$chrom.${my_snps}_unique.rsid \
$temp/combined_uniq.map \
$temp/chr$chrom.${my_snps}_clustering/paintor_chr$chrom.${my_snps}.${my_int}.mcmc
cd $temp
done
done

#Remove empty files (with no results):
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $analysis
cd chr$chrom.${my_snps}_clustering
for b in jam* finemap* paintor*
do
no_lines=$(wc -l $b | cut -d" " -f1 )
if [ $no_lines == 1 ]
then
	rm -rf $b
fi
done
done

#Plot heatmap showing similarity of results.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $temp/chr$chrom.${my_snps}_clustering
echo chr$chrom.${my_snps}_clustering
Rscript $scripts/plot_finemapping_comparison.R ${my_snps}_finemapping_comparison 15 30
done

#Prepare output for integration with other types of annotation.
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320 rs41293864
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $temp/chr$chrom.${my_snps}_clustering
Rscript $scripts/melting.R finemap_${my_snps}_finemapping_comparison.txt finemap_${my_snps}.input
Rscript $scripts/melting.R paintor_${my_snps}_finemapping_comparison.txt paintor_${my_snps}.input
Rscript $scripts/melting.R jam_${my_snps}_finemapping_comparison.txt jam_${my_snps}.input
done


#Due to too high number of NAs in those SNP results, cannot obtain a heatmap for those SNPs.
#Try substituting NAs with 0s in that case (which is incorrect).
for my_snps in rs2218565 rs7146581 rs6473227 rs1249910 rs4713555 rs12188917
do
chrom=$(grep $my_snps $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
cd $temp/chr$chrom.${my_snps}_clustering
rm -rf *.pdf
Rscript $scripts/plot_finemapping_comparison_nona.R ${my_snps}_finemapping_comparison.pdf 15 30
done
