#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/open_targets/
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/open_targets
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015


cd $analysis

#Preparation of Paternoster 2015 GWAS for input. Only lead SNPs.
tail -n +2 $gwas/results.euro.pval.1k | (echo -e "chromosome\tbase_pair_location\tvariant_id\teffect_allele\tother_allele\tbeta\tstandard_error\tp-value" && awk -v OFS="\t" '{print $2, $3, $12, $4, $5, $8, $9, $11}') >paternoster2015_open_targets.tsv

echo -e "chromosome\tbase_pair_location\tvariant_id\teffect_allele\tother_allele\tbeta\tstandard_error\tp-value" >paternoster2015_open_targets_select.tsv 
while read line
do
set -- $line
my_rsid=$1
grep -w $my_rsid paternoster2015_open_targets.tsv >>paternoster2015_open_targets_select.tsv 
done < $gwas/paternoster_2015_index_snps_sorted.txt


#Run POST GAP. Remember to use genome coordinates from GRCh38, although I don't think it matters - look-ups done by rsID.
python POSTGAP.py --database_dir databases --summary_stats paternoster2015_open_targets_select.tsv --disease eczema --output paternoster2015_open_targets_select_annotation.tsv --bayesian

#Within the STOP GAP VM, need to rename databases_dir to databases.
#Currently, the test example provided by authors runs fine but our GWAS data stalls - says that VEP did not respond for variant annotation. Possibly unrecognized rsid used? Could it be that rs1473281 SNP is not recognized by ENSEMBL, however it is not one of our searched for SNPs? Better write to authors then.
