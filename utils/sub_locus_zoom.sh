#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

lzoom=/panfs/panasas01/sscm/qh18484/bin/locuszoom/bin
module add languages/python-2.7.5

cd $PBS_O_WORKDIR

${lzoom}/locuszoom --metal $input_file --markercol $snp_id --pvalcol $pval --refsnp $my_ref --flank $my_flank --pop EUR --build hg19 --source 1000G_March2012 --prefix $my_prefix --no-date --plotonly --cache None