#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/pops
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/pops
gwas=$HOME/scratch/new_gwas/gwas/raw
pops=$HOME/scratch/pops

gwas_name="eczema21_discovery"

cd $analysis

#Activate conda environment with required Python packages
source activate myenv

#Preparation of input GWAS file. 
#The summary statistics file must contain a column containing SNP IDs, p-values, 
#the sample size used per SNP can be included in the analysis and be labeled SNP, P, and N.

tail -n +2 $gwas/results.${gwas_name}.txt | (echo -e "SNP\tP\tN" && awk -v OFS="\t" '{print $12, $11, $6}') >${gwas_name}_pops.tsv

#Download SNP reference for MAGMA run.

#Run MAGMA required by PoPs
cd $pops
$HOME/scratch/magma/magma \
	--bfile 1000G.EUR \
	--gene-annot magma_0kb.genes.annot \
	--pval $analysis/${gwas_name}_pops.tsv ncol=N \
	--gene-model snp-wise=mean \
	--out $gwas_name

sbatch --export=ALL,gwas=$gwas_name $scripts/sub_run_pops.sh

#Annotate PoPs results (gathered into a single folder) with gene names and print into one file.
curl -o "hgnc_complete_set.txt" http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

python $scripts/annotate_pops.py $pops/pops_results/ >${gwas_name}_pops_annotated.txt