#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/pops
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/pops
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
pops=/newhome/qh18484/bin/pops

cd $analysis

#Preparation of input GWAS file. 
#The summary statistics file must contain a column containing SNP IDs, p-values, the sample size used per SNP can be included in the analysis and be labeled SNP, P, and N.

tail -n +2 $gwas/results.euro.pval.1k.dbsnp | (echo -e "SNP\tP\tN" && awk -v OFS="\t" '{print $12, $11, $6}') >paternoster2015_pops.tsv

#Download SNP reference for MAGMA run.

#Run MAGMA required by PoPs
module add apps/magma-1.06
cd $pops
magma \
	--bfile 1000G.EUR \
	--gene-annot magma_0kb.genes.annot \
	--pval $analysis/paternoster2015_pops.tsv ncol=N \
	--gene-model snp-wise=mean \
	--out paternoster2015

qsub -v gwas=paternoster2015 $scripts/sub_run_pops.sh

#Annotate PoPs results (gathered into a single folder) with gene names and print into one file.
curl -o "hgnc_complete_set.txt" ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt

python $scripts/annotate_pops.py $pops/pops_results/ >paternoster2015_pops_annotated.txt