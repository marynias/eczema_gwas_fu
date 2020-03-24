#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -V

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc

eqtlgen=/panfs/panasas01/sscm/qh18484/data/eqtl/eQTLgen
eqtl_f=eQTLsFDR-ProbeLevel.txt.gz

cd $PBS_O_WORKDIR

python $scripts/generate_eqtl_genes_annotation.py $eqtl_results $threshold