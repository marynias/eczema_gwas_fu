#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

module del languages/R-3.4.1-ATLAS
module del languages/python-anaconda3-5.2.0
module add languages/R-3.5.1-ATLAS-gcc-6.1

HOME=/panfs/panasas01/sscm/qh18484
scripts=$HOME/bin/eczema_gwas_fu/bayesian_fm/ref_panel

cd $PBS_O_WORKDIR

Rscript --vanilla $scripts/bigld_gpart_plotld.R ${input_file}.info ${input_file}.geno $my_results $my_chrom $snp_id $my_position $my_interval $my_dataset