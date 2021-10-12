#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/loci_definition
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/loci_definition
genome=$HOME/scratch/new_gwas/genome
gwas=$HOME/scratch/new_gwas/gwas/raw

gwas_name="eczema21_discovery"

cd $analysis

##Gencode GRCH37 gene annotation
head $genome/gencode.v19.annotation_filtered_transcript_appris_final.txt
#Indexed for Giggle
$genome/gencode.v19.annotation_filtered_transcript_appris_final.txt_index
##Gencode GRCH38 gene annotation
head $genome/gencode.v31.annotation_filtered_transcript_appris_final.txt
#Indexed for Giggle
$genome/gencode.v31.annotation_filtered_transcript_appris_final.txt_index

##Link them here
ln -s $genome/gencode.v19.annotation_filtered_transcript_appris_final.txt.bed.gz \
gencode.v19.annotation_filtered_transcript_appris_final.txt.bed.gz 
ln -s $genome/gencode.v19.annotation_filtered_transcript_appris_final.txt_index \
gencode.v19.annotation_filtered_transcript_appris_final.txt_index

ln -s $genome/gencode.v31.annotation_filtered_transcript_appris_final.txt.bed.gz \
gencode.v31.annotation_filtered_transcript_appris_final.txt.bed.gz 
ln -s $genome/gencode.v31.annotation_filtered_transcript_appris_final.txt_index \
gencode.v31.annotation_filtered_transcript_appris_final.txt_index

#Run Giggle to find all the genes overlapping 1 Mbp interval around index SNPs.
cp $gwas/leadSNPs.${gwas_name}_sorted_1Mbp.bed ./
bgzip leadSNPs.${gwas_name}_sorted_1Mbp.bed
temp_genome_index=gencode.v19.annotation_filtered_transcript_appris_final.txt_index
query=leadSNPs.${gwas_name}_sorted_1Mbp.bed.gz
query_results=${gwas_name}_sorted_1Mbp_genes
giggle search -i $temp_genome_index -q $query -v -o  > $query_results
python $HOME/bin/eczema_gwas_fu/annotation/giggle/process_giggle_output.py \
$query_results ${query_results}_processed

#Run the following script to annotate RSID with cytobands.
bedtools intersect -loj -a $gwas/leadSNPs.${gwas_name}.bed -b $genome/ucsc_hg19_cytoBand.bed > $gwas/leadSNPs.${gwas_name}_cytoBand.bed

