#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015

cd $PBS_O_WORKDIR

python $scripts/gene_generate_eqtl_input_pvals_zscore_by_rsid.py --tab QTLsFDR-ProbeLevel_genes_eczema3Mbp.txt --gene 17 \
--proc $gwas/paternoster_2015_index_snps_sorted.txt --interval $my_int --snp $my_snps \
--genes $HOME/analysis/colocalization/coloc/sun_pqtl/${my_snps}_${my_int}.gene_names_abbrv --chrom 3 --pos 4 --pval 1 --ident 2 --zscore 11