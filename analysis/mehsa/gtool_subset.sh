#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -V

########Subset Geno and sample files to the list of specified individuals

# on compute node, change directory to 'submission directory':
cd $analysis_dir

var_data=/panfs/panasas01/sscm/qh18484/twinsuk_var/released/2017-08-16/data
analysis_dir=/panfs/panasas01/sscm/qh18484/analysis/mehsa

for i in {01..22}
do
gtool -S --g $var_data/gen/data.chr${i}.gen \
--s $var_data/sample/data.chr${i}.sample \
--og $analysis_dir/data_subset.chr${i}.gen \
--os $analysis_dir/data_subset.chr${i}.sample \
--sample_id $analysis_dir/examples/Sample_ID.txt \
--log $analysis_dir/subset_all_chr${i}.log 
done