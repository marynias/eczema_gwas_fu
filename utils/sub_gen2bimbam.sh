#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

#Convert from GEN to BIMBAM

# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR

cat $gen_file | awk -v s=$number_of_samples \
'{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) \
printf "," $(i*3+3)*2+$(i*3+4); \
printf "\n" }' > $bb_file