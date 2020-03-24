#!/bin/bash
dataset=/panfs/panasas01/sscm/qh18484/working/data/Datasets/RefPanel/UK10k
analysis=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/UK10k
input_vcf=UK10K_COHORT.20160215.sites.vcf
cd $dataset

#Download the VCF file with UK10k dataset, its md5 file and README.
wget ftp://ngs.sanger.ac.uk/production/uk10k/UK10K_COHORT/REL-2012-06-02/UK10K_COHORT.20160215.sites.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/uk10k/UK10K_COHORT/REL-2012-06-02/UK10K_COHORT.20160215.sites.vcf.gz.md5
wget ftp://ngs.sanger.ac.uk/production/uk10k/UK10K_COHORT/REL-2012-06-02/README
#Generate md5sum hash to see if file not corrupted in transfer
md5sum -c ${input_vcf}.gz.md5
#File is OK.

cd $analysis
#Extract the file to working folder on the HPC.
gunzip -c $dataset/${input_vcf}.gz > $analysis/${input_vcf}

#Create a fake panel file with which to use the CalcLD_1KG_VCF.py script from PAINTOR suite.
grep -v "##" $input_vcf | head -1 | cut -f10- | tr '\t' '\n' | awk '{print $1 "\t" "EUR" "\t" "EUR" "\t" "NA"}' | sed '1s/^/sample\tpop\tsuper_pop\tgender\n/'

