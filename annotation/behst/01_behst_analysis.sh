#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/behst
scripts=$HOME/bin/eczema_gwas_fu/annotation/behst

#Identify the GO enrichment in genes whose regulatory regions interact with regions around focal SNPs. 
#Use 1Mbp and 3Mbp, 100kbp and 250kbp around focal SNP - files created during Sun 2018 pQTL colocalisation.
#Generate bed files based on bigld intervals.
cd $analysis

#Concantenate all regions for a given interval to be used in BEHST analysis.
cat *_100kbp.bed >all_100kbp.bed
cat *_250kbp.bed >all_250kbp.bed
cat *_1Mbp.bed >all_1Mbp.bed
cat *_3Mbp.bed >all_3Mbp.bed
cat *bigld.bed >all_bigld.bed

#Sort the files
for a in all*bed
do
sort -g -k1,2 $a >temp
mv temp $a
done

#The bed files contain redundant regions: secondary effects and replicated loci. Keep only one representative interval per locus in each file.
for a in all*bed
do
Rscript $scripts/merge_redundant_regions.R $a paternoster2015_representative_snps.txt ${a%.bed}_nonred.bed
done

#add chr to column 1 to match the BEHST requirements.
for a in all*bed
do
sed -i s/^/chr/ $a
done

#Run behst
for a in all*_nonred.bed
do
behst $a -d /panfs/panasas01/sscm/qh18484/bin/behst
done