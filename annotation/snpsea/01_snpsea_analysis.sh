##!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
snpsea=$HOME/bin/snpsea
scripts=$HOME/bin/eczema_gwas_fu/annotation/snpsea
analysis=$HOME/analysis/annotation/snpsea
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
suffix1=paternoster2015_index_snps
cd $analysis
#Prepare SNPSea input - both only index SNPs and SNPs within interval.
#Interval
(printf "CHR\tPOS\tSNP\tP\n"; cat $data_manipulation/interval_r2_0.2_1k.gwas4d | cut -f1,2,3,6 | sed 's/^/chr/' | sort -k1,1V -k2,2n ) >interval_r2_0.2_1k.snpsea
#Index SNP only.
cut -f1 $gwas/paternoster_2015_index_snps_sorted.txt >paternoster2015_index_snps
(printf "CHR\tPOS\tSNP\tP\n"; grep -w -F -f  paternoster2015_index_snps interval_r2_0.2_1k.snpsea) >paternoster2015_index_snps.snpsea

$snpsea/bin/snpsea --snps $snpsea/Red_blood_cell_count-Harst2012-45_SNPs.gwas \
--gene-matrix  $snpsea/GeneAtlas2004.gct.gz \
--gene-intervals  $snpsea/NCBIgenes2013.bed.gz \
--snp-intervals  $snpsea/TGP2011.bed.gz \
--null-snps $snpsea/Lango2010.txt.gz \
--out out \
--slop 10e3 \
--threads 4 \
--null-snpsets 1000 \
--min-observations 100 \
--max-iterations 1e7 


qsub -v my_gwas=$snpsea/Red_blood_cell_count-Harst2012-45_SNPs.gwas,gene_matrix=$snpsea/GeneAtlas2004.gct.gz,output=out $scripts/sub_snpsea.sh
qsub -v input=out $scripts/sub_snpsea_viz.sh

qsub -v my_gwas=paternoster2015_index_snps.snpsea,gene_matrix=$snpsea/GeneAtlas2004.gct.gz,output=${suffix1}_GeneAtlas2004 $scripts/sub_snpsea.sh
qsub -v my_gwas=paternoster2015_index_snps.snpsea,gene_matrix=$snpsea/ImmGen2012.gct.gz,output=${suffix1}_ImmGen2012 $scripts/sub_snpsea.sh
qsub -v my_gwas=paternoster2015_index_snps.snpsea,gene_matrix=$snpsea/FANTOM2014.gct.gz,output=${suffix1}_FANTOM2014 $scripts/sub_snpsea.sh
qsub -v my_gwas=paternoster2015_index_snps.snpsea,gene_matrix=$snpsea/GO2013.gct.gz,output=${suffix1}_GO2013 $scripts/sub_snpsea.sh

#Sort the results
for a in ${suffix1}_GeneAtlas2004 ${suffix1}_ImmGen2012 ${suffix1}_FANTOM2014 ${suffix1}_GO2013
do
(head -1 $a/condition_pvalues.txt; tail -n +2 $a/condition_pvalues.txt | sort -t$'\t' -k2n) > $a/condition_pvalues_sorted.txt 
done