#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/loci_definition/eqtl_catalogue
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/loci_definition
cytoband=$HOME/analysis/network/AD

cd $analysis

#Note: can use MendelVar to generate locus intervals, however if SNP is missing in the database,
#will have to find its proxy manually. May want to use the following synonym list:
$HOME/analysis/annotation/data_manipulation/rsid_synonyms.txt

##Reference files from MendelVar analysis
##Gencode GRCH37 gene annotation
head $HOME/MendelVar/genome/gencode.v19.annotation_filtered_transcript_appris_final.txt
#Indexed for Giggle
$HOME/MendelVar/genome/gencode.v19.annotation_filtered_transcript_appris_final.txt_index
##Gencode GRCH38 gene annotation
head $HOME/MendelVar/genome/ggencode.v31.annotation_filtered_transcript_appris_final.txt
#Indexed for Giggle
$HOME/MendelVar/genome/gencode.v31.annotation_filtered_transcript_appris_final.txt_index

##Link them here
ln -s $HOME/MendelVar/genome/gencode.v19.annotation_filtered_transcript_appris_final.txt.bed.gz \
gencode.v19.annotation_filtered_transcript_appris_final.txt.bed.gz 
ln -s $HOME/MendelVar/genome/gencode.v19.annotation_filtered_transcript_appris_final.txt_index \
gencode.v19.annotation_filtered_transcript_appris_final.txt_index

ln -s $HOME/MendelVar/genome/gencode.v31.annotation_filtered_transcript_appris_final.txt.bed.gz \
gencode.v19.annotation_filtered_transcript_appris_final.txt.bed.gz 
ln -s $HOME/MendelVar/genome/gencode.v31.annotation_filtered_transcript_appris_final.txt_index \
gencode.v19.annotation_filtered_transcript_appris_final.txt_index
#Run Giggle to find all the genes overlapping 3 Mbp interval around index SNPs.
bgzip paternoster_2015_index_snps_sorted_1Mbp.bed
temp_genome_index=gencode.v19.annotation_filtered_transcript_appris_final.txt_index
query=paternoster_2015_index_snps_sorted_1Mbp.bed.gz
query_results=paternoster_2015_index_snps_sorted_1Mbp_genes
giggle search -i $temp_genome_index -q $query -v -o  > $query_results
python $HOME/bin/eczema_gwas_fu/annotation/giggle/process_giggle_output.py \
paternoster_2015_index_snps_sorted_1Mbp_genes latest

#Run the following script to annotate RSID with cytobands.
intersectBed -loj -a paternoster_2015_index_snps.bed -b $cytoband/ucsc_hg19_cytoBand.bed >paternoster_2015_index_snps_cytoband.bed

