#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

########Subset Geno and sample files to the list of specified individuals
# on compute node, change directory to 'submission directory':
cd $PBS_O_WORKDIR

gtool -S --g $InGen \
--s $InSample \
--og $OutGen \
--os $OutSample \
--sample_id $Subset \
--inclusion $To_include \
--log $LogFile \