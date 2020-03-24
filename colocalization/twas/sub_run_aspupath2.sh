#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas
twas=$HOME/bin/fusion_twas-master
my_gwas=results.euro.pval.1k

cd $PBS_O_WORKDIR

Rscript aSPUpath2.R \
--sumstats ./Example/example.stat.rds \
--out ./Example/example_res.rds \
--weights ./WEIGHTS/CMC.BRAIN.RNASEQ.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld ./LDREF/1000G.EUR. \
--pathway_list ./Example/example_GOCC.txt

Rscript $twas/FUSION.post_process.R \
--sumstats $analysis/${my_gwas%.1k}.chr${chrom} \
--input $analysis/${tissue}_${my_gwas%.1k}.twas.chr${chrom}.top \
--out $analysis/${tissue}_${my_gwas%.1k}.twas.chr${chrom}.analysis \
--ref_ld_chr $twas/LDREF/1000G.EUR. \
--chr $chrom \
--plot --plot_corr --plot_eqtl --plot_individual --plot_scatter --report --save_loci --locus_win 100000