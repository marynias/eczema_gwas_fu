#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

bgenix=/panfs/panasas01/sscm/qh18484/bin/bgenix/bin
ukbiobank=/panfs/panasas01/sscm/qh18484/data/ukbiobank

cd $PBS_O_WORKDIR

%%	rsid=do_wymiany
%%	chrom=do_wymiany
%%	pos=do_wymiany
%%	interval=do_wymiany

start="$((pos - 1500000))"
end="$((pos + 1500000))"

printf -v chrom_2 "%02d" $chrom

$bgenix/bgenix -g $ukbiobank/dosage_bgen/data.chr${chrom_2}.bgen -incl-range ${chrom_2}:${start}-${end} -vcf | gzip -c > chr1${chrom}_${rsid}_${interval}_ukbb.vcf.gz
