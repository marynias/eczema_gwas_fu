#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00

#Convert from GEN to BIMBAM

# on compute node, change directory to 'submission directory':
analysis_dir=/panfs/panasas01/sscm/qh18484/analysis/mehsa
cd $analysis_dir

number_of_samples=672
for i in {01..22}
do
gen_file=$analysis_dir/data_subset.chr${i}.gen
bb_file=$analysis_dir/chr${i}_gen_test2.bimbam 
cat $gen_file | awk -v s=$number_of_samples \
'{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) \
printf "," $(i*3+3)*2+$(i*3+4); \
printf "\n" }' > $bb_file
done