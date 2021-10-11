#!/bin/bash
set -e -u -o pipefail
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/new_gwas/colocalisation/eqtl_catalogue
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue

cd $analysis

###!! Use R 4.0 for coloc (normally using R 3.6)

####!!!! Remeber to update the s parameter for eczema - proportion of cases.
####eQTL catalogue requires GRCh38 intervals, while eQTLgen and Sun require GRCh37.
####However, merging by rsids.


##Loop over all lead RSIDs and run colocalisation analysis for transcripts:
for my_r in $(cat $my_bed_file | cut -f4)
do
	echo $my_r
	#Rscript --vanilla $scripts/eQTL_Catalogue_coloc_transcript.R $my_r
	qsub -v my_rsid=$my_r $scripts/sub_run_eqtl_transcript.sh
done


#In the end, pulling data from the API is too slow due to 
#frequent stalling - only 20 colocalisations
#done over 25 minutes. Going to download gene eQTLs and run locally.

#Download eQTL Catalog files (ge) to disk.
cd catalog
sh $scripts/get_eQTL_catalog_files.sh

#Run coloc using files on disk.
for my_r in $(cat $my_bed_file | cut -f4)
do
	echo $my_r
	#Rscript --vanilla $scripts/eQTL_Catalogue_coloc_local.R $my_r
	qsub -v my_rsid=$my_r $scripts/sub_run_eqtl_local.sh
done


#Run coloc using Sun et al. (2018) and eQTLgen from MRC IEU GWAS db.
#All loci in single job, as running quick.
#Run it on my local computer.
#Make sure that all the stop finish running.

qsub -v $scripts/sub_run_ieugwasr.sh

#Merge all the results from different loci into one table for gxp and tx.
#Remove empty files with no rows (just header) before join.
Rscript --vanilla $scripts/join_all_tables.R