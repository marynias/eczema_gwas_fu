#!/bin/bash
HOME=/mnt/storage/home/qh18484
analysis=$HOME/scratch/new_gwas/colocalisation/eqtl_catalogue
scripts=$HOME/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue
loci=$HOME/scratch/new_gwas/loci_definition
utils=$HOME/bin/eczema_gwas_fu/utils
gwas=$HOME/scratch/new_gwas/gwas/raw

cd $analysis

gwas_name="eczema21_discovery"

#mkdir -p GWAS_intervals

####!!!! Remeber to update the s parameter for eczema - proportion of cases.
####eQTL catalogue requires GRCh38 intervals, while eQTLgen and Sun require GRCh37.
####However, merging by rsids.

#Pulling data from eQTL Catalogue API is too slow due to 
#frequent stalling - only 20 colocalisations
#done over 25 minutes. Going to download eQTLs and run all locally.

#Download eQTL Catalog files (ge) to disk.
#sh $scripts/get_eQTL_catalog_files.sh

#Create a set of files holding GWAS Summary stats in the selected GWAS intervals.
while read line
do
echo $line >temp
rsid=$(cat temp | cut -d" " -f4)
python $utils/subset_interval_chr.py --tab $gwas/results.${gwas_name}.txt --header_tab Y --bed temp --pos_tab 3 --chr_tab 2 --output ${rsid}_${gwas_name}.gwas
mv ${rsid}_${gwas_name}.gwas GWAS_intervals/
done < $loci/${gwas_name}_interval.bed

my_bed_file=$loci/${gwas_name}_interval.bed
#Run coloc using files on disk.
for my_r in $(cat $my_bed_file | cut -f4)
do
	#Rscript --vanilla $scripts/eQTL_Catalogue_coloc_local.R $my_r $loci/${gwas_name}_interval_mapped.bed \
	#$loci/${gwas_name}_sorted_1Mbp_genes_processed $gwas_name
	sbatch --export=ALL,my_rsid=$my_r,hg38=$loci/${gwas_name}_interval_mapped.bed,gene_list=$loci/${gwas_name}_sorted_1Mbp_genes_processed,gwas=$gwas_name $scripts/sub_run_eqtl_local.sh 
done

##### Rewrite to use local files from BC4, rather than calls to API.
#Run coloc using Sun et al. (2018) and eQTLgen from MRC IEU GWAS db.
#All loci in single job, as running quick.
#Run it on my local computer.
sbatch --export=ALL,my_interval=$loci/${gwas_name}_interval.bed,gene_list=$loci/${gwas_name}_sorted_1Mbp_genes_processed,gwas=$gwas_name $scripts/sub_run_ieugwasr.sh

#Merge all the results from different loci into one table for gxp and tx.
#Remove empty files with no rows (just header) before join.
Rscript --vanilla $scripts/join_all_tables.R