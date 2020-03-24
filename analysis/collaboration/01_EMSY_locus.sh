#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/collaboration
scripts=$HOME/bin/eczema_gwas_fu/annotation/collaboration
skin_data=$HOME/data/eqtl/twinsuk/skin
var_data=$HOME/data/eqtl/twinsuk/gen
sample_data=$HOME/data/eqtl/twinsuk/sample
#Obtain RPKM gene expression in the skin, and genotypes for gene "ENSG00000158636" (EMSY) and SNP rs7101927 in TwinsUK.
cd $analysis
head -1 $skin_data/data.rpkm >rs7101927_EMSY.rpkm
grep -w "ENSG00000158636" $skin_data/data.rpkm >>rs7101927_EMSY.rpkm
grep -w "rs7101927" $var_data/data.chr11.gen >rs7101927.gen
cp $sample_data/data.chr11.sample rs7101927.sample
plink --gen rs7101927.gen --sample rs7101927.sample --recodeAD --allow-extra-chr --out rs7101927