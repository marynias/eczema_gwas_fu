#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR

%%	input_file=do_wymiany
%%  optional=" "
plink --vcf $input_file $optional --out ${input_file%.vcf*} > ${input_file%.vcf*}.log