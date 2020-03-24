#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/KGGSeq
scripts=$HOME/bin/eczema_gwas_fu/annotation/KGGSeq
prog_dir=$HOME/bin/kggseq10hg19
data_manipulation=$HOME/analysis/annotation/data_manipulation

cd $analysis

for a in $data_manipulation/r2_0.2_1k.gwas4d $data_manipulation/paternoster_2015_index_snps_sorted_1Mbp.gwas4d \
$data_manipulation/paternoster_2015_index_snps_sorted_3Mbp.gwas4d
do
cat $a | awk -v OFS="\t" '{print $1, $2, $2, $4, $5, $3}' >${a%.gwas4d}.kggseq
java -jar $prog_dir/kggseq.jar --annovar-file ${a%.gwas4d}.kggseq --db-score dbncfp_known --regulatory-causing-predict all
gunzip kggseq1.flt.txt.gz
mv kggseq1.flt.txt ${a%.gwas4d}.kggseq_out 
done 


