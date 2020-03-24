#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=12:00:00

paintor=/panfs/panasas01/sscm/qh18484/bin/PAINTOR_V3.0/PAINTOR_Utilities

cd $PBS_O_WORKDIR

%%	my_file=do_wymiany
%%	my_reference=do_wymiany
%%	my_map=do_wymiany
%%	chrom=do_wymiany
%%	my_eff=do_wymiany
%%	my_alt=do_wymiany
%%	my_pop=do_wymiany
%%	my_zscore=do_wymiany
%%	my_pos=do_wymiany

python $paintor/CalcLD_1KG_VCF_keep_ambiguous.py \
--locus $my_file \
--reference $my_reference \
--map $my_map \
--effect_allele $my_eff \
--alt_allele $my_alt \
--population $my_pop \
--Zhead $my_zscore \
--out_name ${my_file%.locus} \
--position $my_pos
