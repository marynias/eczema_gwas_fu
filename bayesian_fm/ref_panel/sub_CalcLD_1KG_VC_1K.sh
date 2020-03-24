#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00

paintor=/panfs/panasas01/sscm/qh18484/bin/PAINTOR_V3.0/PAINTOR_Utilities

cd $PBS_O_WORKDIR

%%	my_file=do_wymiany
%%	chrom=do_wymiany

python $paintor/CalcLD_1KG_VCF_keep_ambiguous.py \
--locus $my_file \
--reference ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz \
--map integrated_call_samples_v3.20130502.EUR.panel \
--effect_allele A1 \
--alt_allele A0 \
--population EUR \
--Zhead Zscore \
--out_name ${my_file%.locus*} \
--position pos
