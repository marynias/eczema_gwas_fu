#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
snpsea=$HOME/bin/snpsea

cd $PBS_O_WORKDIR

$snpsea/bin/snpsea --snps $my_gwas \
--gene-matrix  $gene_matrix \
--gene-intervals  $snpsea/NCBIgenes2013.bed.gz \
--snp-intervals  $snpsea/TGP2011.bed.gz \
--null-snps $snpsea/Lango2010.txt.gz \
--out $output \
--slop 10e3 \
--threads 2 \
--null-snpsets 0 \
--min-observations 100 \
--max-iterations 1e7 