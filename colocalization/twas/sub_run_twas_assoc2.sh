#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas/gtex
twas=$HOME/bin/fusion_twas-master
my_gwas=results.euro.pval.1k

cd $PBS_O_WORKDIR

Rscript $twas/FUSION.assoc_test.R \
--sumstats $analysis/${my_gwas%.1k}.chr${chrom} \
--weights ${input}/${tissue}.P01.pos \
--weights_dir ${input} \
--ref_ld_chr $twas/LDREF/1000G.EUR. \
--chr $chrom \
--perm 100 \
--out ${input}/${tissue}_out/${my_gwas%.1k}.twas.chr${chrom}
