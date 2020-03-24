#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/data_manipulation
scripts=$HOME/bin/eczema_gwas_fu/annotation/data_manipulation
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
utils=$HOME/bin/eczema_gwas_fu/utils
###Annotation and subsetting of datasets used in GWAS fine-mapping analysis.
cd $analysis

#Add OR ratio to GWAS results.
cd $gwas
for a in results.euro.pval.1k results.euro.pval.all.1k results.euro.pval.all.ukbiobank results.euro.pval.ukbiobank
do
Rscript $scripts/add_log_odds.R $a ${a}.or
done 

#Add chromosome info (make the file into bed file) to the file with r2-based (r2>=0.2) interval generated using 1k EUR data.
while read line
do
set -- $line
my_start=$1
my_end=$2
my_rsid=$3
my_chrom=$(grep -w $my_rsid $gwas/paternoster_2015_index_snps_sorted.txt | cut -f2)
echo -e $my_chrom'\t'$my_start'\t'$my_end'\t'$my_rsid >>interval_r2_0.2_1k.bed
done < interval_r2_0.2_1k.txt 

#Generate the same for fixed 1Mbp and 3Mbp intervals (copied over from the Sun2018 pQTL analysis).
#Also, create "chr"-padded versions for chromosome names.
for a in *.bed
do
sed 's/^/chr/' $a >${a%.bed}_padded.bed
done

#Create non-redundant datasets which does not include replicated loci and secondary signals.
for a in *bed
do
python $utils/subset_snps.py --tab $a --snp_list paternoster2015_non_redunant_list.txt --header_tab N --pos_tab 4 --output ${a%.bed}_non_redundant.bed
done

#Create VCF-style input for GWAS4D
for a in interval_r2_0.2_1k.bed interval_r2_0.2_1k_non_redundant.bed paternoster_2015_index_snps_sorted_1Mbp.bed paternoster_2015_index_snps_sorted_1Mbp_non_redundant.bed paternoster_2015_index_snps_sorted_3Mbp.bed paternoster_2015_index_snps_sorted_3Mbp_non_redundant.bed
do
python $utils/subset_interval_chr.py --bed $a --tab $gwas/results.euro.pval.1k --header_tab Y --pos_tab 3 --chr_tab 2 \
 --output ${a%.bed}.gwas
done 

#Subset BED-style interval with the second to last column being SNP ID, and the last column being index SNP.
python $utils/subset_interval_chr_map.py --bed interval_r2_0.2_1k_nodups.bed --tab $gwas/results.euro.pval.1k --header_tab Y --pos_tab 3 --chr_tab 2 \
 --output interval_r2_0.2_1k_nodups.gwas

tail -n +2 interval_r2_0.2_1k_nodups.gwas | awk -v OFS="\t" '{print $2, $3, $3, $12, $13}' >interval_r2_0.2_1k_nodups.map

for a in *gwas
do
cat $a | awk -v OFS="\t" '{print $2, $3, $12, $4, $5, $11}' | tail -n +2 >${a}4d
done

#Create a list of SNPs for GoDMC.
for a in interval_r2_0.2_1k.gwas paternoster_2015_index_snps_sorted_3Mbp.gwas paternoster_2015_index_snps_sorted_1Mbp.gwas
do
cut -f1,12 $a | sed 's/^/chr/' > ${a%.gwas}.godmc
done

#Bed files containing all individual location of SNPs in the analysis.
for a in interval_r2_0.2_1k.gwas4d paternoster_2015_index_snps_sorted_3Mbp.gwas4d paternoster_2015_index_snps_sorted_1Mbp.gwas4d
do
cat $a | awk -v OFS="\t" '{print $1, $2, $2, $3}' > ${a%.gwas4d}_individual.bed
done

#Generate list of all the SNPs in the interval.
cut -f4 interval_r2_0.2_1k_individual.bed >interval_r2_0.2_1k.snps


#Generate VCF input for FATHMM-XF and DeepSEA. Need to generate two versions: with effect allele as the referenec allele (ref file) and with non-effect allele as the reference allele (alt file)
cat $gwas/results.euro.pval.1k | awk -v OFS="\t" '{print $2, $3, $12, $4, $5}' >results.euro.pval.1k.ref.vcf
cat $gwas/results.euro.pval.1k | awk -v OFS="\t" '{print $2, $3, $12, $5, $4}' >results.euro.pval.1k.alt.vcf

#Subset to SNPs of interest.
for a in ref alt 
do
python $utils/subset_interval_chr.py --bed interval_r2_0.2_1k.bed --tab results.euro.pval.1k.${a}.vcf --header_tab N --pos_tab 2 --chr_tab 1 \
 --output interval_r2_0.2_1k_${a}.vcf
done

#Convert BED to granges format.
for a in interval_r2_0.2_1k.bed interval_r2_0.2_1k_individual.bed 
do
cat $a | sed 's/^/chr/' | cut -f1-3 | awk '{print $1":"$2"-"$3}' >${a%bed}granges
done

#Generate a ref file listing gene name, chromosome and transcript start and end - based on Gencode
python $scripts/gene_start_finish.py ./gtf/gencode.v19.annotation.gtf >gencode.v19.annotation.se
sed -i 's/"//g' gencode.v19.annotation.se

#Generate a ref file listing gene's ensembl ID; name, chromosome and transcript start and end - based on Gencode.
python $scripts/gene_start_finish_ensembl.py ./gtf/gencode.v19.annotation.gtf >gencode.v19.annotation.ensembl
sed -i 's/"//g' gencode.v19.annotation.ensembl

#Generate a ref file listing gene name, chromosome and transcript start and end - based on Ensembl.
#Download the annotation gtf file from ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
python $scripts/gene_start_finish_ensembl.py ./gtf/Homo_sapiens.GRCh37.87.gtf >Homo_sapiens.GRCh37.87.ensembl
sed -i 's/"//g' Homo_sapiens.GRCh37.87.ensembl

#Filter the HUGO gene name reference file to remove withdrawn entries:
grep -v "entry withdrawn" hugo_synonyms_ids2.txt | grep -v "symbol withdrawn" | sed 's/"//g' >hugo_synonyms_ids2_filtered.txt

#Create a JSON-saved dictionary file listing all possible synonyms of a given Ensembl ENSG ID:
#key - EnsemblID, value-synonym

#Create a JSON-saved dictionary file listing all possible synonyms of a HUGO gene name:
#key - HUGO gene name, value-synonym
#Remove Non-ASCII characters.
perl -pi -e 's/[^[:ascii:]]//g' hugo_synonyms_ids2_filtered.txt
python3 $scripts/gene_json_dictionary.py hugo_synonyms_ids2_filtered.txt
#Remove question marks left over by non-ASCII character substitution.
sed -i 's/?-//g' hugo_synonyms_ids2_filtered.ensembl
sed -i 's/?//g' hugo_synonyms_ids2_filtered.ensembl
sed -i 's/?-//g' hugo_synonyms_ids2_filtered.hugo
sed -i 's/?//g' hugo_synonyms_ids2_filtered.hugo
sed -i 's/?//g' hugo_synonyms_ids2_filtered.hugo.table
sed -i 's/?//g' hugo_synonyms_ids2_filtered.ensembl.table
#Edit the files by hand to remove the first few entries with just punctuation marks.
#Create version all in lowercase or uppercase.
tr '[:upper:]' '[:lower:]' < hugo_synonyms_ids2_filtered.hugo > hugo_synonyms_ids2_filtered.hugo.lower
tr '[:upper:]' '[:lower:]' < hugo_synonyms_ids2_filtered.ensembl > hugo_synonyms_ids2_filtered.ensembl.lower
tr '[:upper:]' '[:lower:]' < hugo_synonyms_ids2_filtered.ensembl.table > hugo_synonyms_ids2_filtered.ensembl.table.lower
tr '[:upper:]' '[:lower:]' < hugo_synonyms_ids2_filtered.hugo.table > hugo_synonyms_ids2_filtered.hugo.table.lower

tr '[:lower:]' '[:upper:]' < hugo_synonyms_ids2_filtered.hugo > hugo_synonyms_ids2_filtered.hugo.upper
tr '[:lower:]' '[:upper:]'  < hugo_synonyms_ids2_filtered.ensembl > hugo_synonyms_ids2_filtered.ensembl.upper
tr '[:lower:]' '[:upper:]' < hugo_synonyms_ids2_filtered.ensembl.table > hugo_synonyms_ids2_filtered.ensembl.table.upper
tr '[:lower:]' '[:upper:]'  < hugo_synonyms_ids2_filtered.hugo.table > hugo_synonyms_ids2_filtered.hugo.table.upper


#Create a dictionary with keys = dbsnp current ref rsid, value = rsid alias
python $scripts/snp_json_dictionary.py $HOME/working/data/dbSNP/SNP-Base_files_merged/aliastable.final.txt

#Create a bim file from the whole dbSNP database
zcat $HOME/working/data/dbSNP/SNP-Base_files_merged/dbsnptable.AFchecked.cleaned.Sept2018.txt.gz | tr ' ' '\t'  |  sed 's/:/\t/g' | awk -v OFS="\t" '{print $4, $1, '0', $5, $8, $9}'  >dbsnptable.AFchecked.cleaned.Sept2018.bim
cat dbsnptable.AFchecked.cleaned.Sept2018.bim | tail -n +1 >temp
mv temp dbsnptable.AFchecked.cleaned.Sept2018.bim

zcat $HOME/working/data/dbSNP/SNP-Base_files_merged/dbsnptable.AFchecked.allconcat.txt.gz | tr ' ' '\t' | awk -v OFS="\t" '{print $6, $1, '0', $7, $4, $5}' >dbsnptable.AFchecked.allconcat.bim
cat dbsnptable.AFchecked.allconcat.bim | tail -n +1 >temp
mv temp dbsnptable.AFchecked.allconcat.bim