#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

HOME=/panfs/panasas01/sscm/qh18484
priority_pruner=$HOME/bin/PriorityPrunner

cd $PBS_O_WORKDIR

%%	tfile=do_wymiany
%%	snp_table=do_wymiany
%%	interval=do_wymiany
%%	r2=do_wymiany
%%	min_maf=do_wymiany
%%	min_snp_call_rate=do_wymiany
%%	output=do_wymiany

java -jar $priority_pruner/PriorityPruner.jar --tfile $tfile --snp_table $snp_table \
--r2 $r2 --max_distance $interval --min_maf $min_maf \
 --min_snp_callrate $min_snp_call_rate --out $output