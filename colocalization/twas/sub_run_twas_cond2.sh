#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas/gtex
twas=$HOME/bin/fusion_twas-master
my_gwas=results.euro.pval.1k

cd $PBS_O_WORKDIR

no_genes=$(wc -l $input/${tissue}.P01.pos | cut -d" " -f1)
cat ${input}/${tissue}_out/${my_gwas%.1k}.twas.chr${chrom} | awk -v genes="$no_genes" 'NR == 1 || $20 < 0.05/genes' >${input}/${tissue}_out/${my_gwas%.1k}.twas.chr${chrom}.top

Rscript $twas/FUSION.post_process.R \
--sumstats $analysis/${my_gwas%.1k}.chr${chrom} \
--input ${input}/${tissue}_out/${my_gwas%.1k}.twas.chr${chrom}.top \
--out ${input}/${tissue}_out/${my_gwas%.1k}.twas.chr${chrom}.analysis \
--ref_ld_chr $twas/LDREF/1000G.EUR. \
--chr $chrom \
--plot --locus_win 100000