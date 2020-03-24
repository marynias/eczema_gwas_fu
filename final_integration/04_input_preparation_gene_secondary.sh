#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/final_integration/gene_secondary
scripts=$HOME/bin/eczema_gwas_fu/final_integration
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
data_manipulation=$HOME/analysis/annotation/data_manipulation
data_manipulation_scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
datasets=$HOME/working/data/Datasets

cd $analysis

#Annotate all the genes-only processed files with locus IDs established in the lookups processed files.
python $scripts/annotate_genes_locus_id.py 
