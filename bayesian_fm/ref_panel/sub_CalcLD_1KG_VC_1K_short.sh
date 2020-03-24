#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00

paintor=/panfs/panasas01/sscm/qh18484/bin/PAINTOR_V3.0/PAINTOR_Utilities
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes

cd $analysis
python $paintor/CalcLD_1KG_VCF_keep_ambiguous.py \
--locus chr1_short.locus \
--reference ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--map integrated_call_samples_v3.20130502.ALL.panel \
--effect_allele A1 \
--alt_allele A0 \
--population EUR \
--Zhead Zscore \
--out_name chr1_ambig_short.ld \
--position pos
