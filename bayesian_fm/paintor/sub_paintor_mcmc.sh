#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -V

# on compute node, change directory to 'submission directory':
set -vx
cd $PBS_O_WORKDIR
PAINTOR=$HOME/bin/PAINTOR_V3.0
#Run PAINTOR without use of any annotation, for bigger number of SNP states considered.

$PAINTOR/PAINTOR -input $input \
-Zhead Zscore  \
-LDname ld \
-in $PWD \
-out $PWD \
-max_causal $max_causal \
-MI 10000 \
-mcmc \
-annotations dummy