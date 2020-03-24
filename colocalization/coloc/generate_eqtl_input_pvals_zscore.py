#! /usr/bin/env python
#Generate eQTL input table for coloc, if we only have p-values available.
#eQTL results can be gzipped.

#! /usr/bin/env python
import os, sys, re, argparse
import gzip
from collections import defaultdict as dd
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with eQTL summary stats')
ap.add_argument('--proces',required=True,type=str,help="File with markers id and location (chrom, pos) for index SNPs to be used in colocalisation analysis")
ap.add_argument('--interval',required=True,type=int,help='Total interval size centred on index SNP [in kbp] ')   
ap.add_argument('--chrom',required=True,type=int,help='Column with chrom information in the --tab file')  
ap.add_argument('--pos',required=True,type=int,help='Column with pos information in the --tab file') 
ap.add_argument('--pval',required=True,type=int,help='Column with pval information in the --tab file')   
ap.add_argument('--ident',required=True,type=int,help='Column with snp id information in the --tab file')  
ap.add_argument('--zscore',required=True,type=int,help='Column with Zscore information in the --tab file')    

args = ap.parse_args()
#Initialize
tab = args.tab
interval = args.interval
proces = args.proces
chrom = args.chrom
pos = args.pos
pval = args.pval
ident = args.ident
zscore = args.zscore

my_loci = {}
files_dict = {}

#Read in loci to be processed and create appropriate output files.
with open (proces) as proces_h:
	for line in proces_h:
		lines = line.strip().split()  
		my_loci[lines[0]] = [int(lines[1]), int(lines[2])]
		filename = lines[0] + "_" + str(interval) + "kbp.eqtl"
		files_dict.setdefault(filename, open(filename,'w')) 
		files_dict[filename].write("rsid" + "\t" + "chrom" + "\t" + "pos" + "\t" + "pval" + "\t" + "zscore" + "\n")

#Read in gene names to be processed

if tab.endswith(".gz"):
	file_stream = gzip.open(tab, 'rt')
else:
	file_stream = open(tab, 'r')

def evaluate_snp(my_id,my_chrom,my_pos,my_pval, my_zscore):
	for snp_index in my_loci.keys():
		snp_chrom = my_loci[snp_index][0]
		if snp_chrom == my_chrom:
			lower_bound = my_loci[snp_index][1] - (float(interval) * 1000 / 2)
			upper_bound = my_loci[snp_index][1] + (float(interval) * 1000 / 2)
			if (my_pos >= lower_bound and my_pos <= upper_bound):
				filename =  snp_index + "_" + str(interval) + "kbp.eqtl"
				files_dict[filename].write(my_id + "\t" + str(my_chrom) + "\t" + str(my_pos) + "\t" + str(my_pval) + "\t" + str(my_zscore) + "\n")
			else:
				continue
		else:
			continue

#Scan eQTL file for target SNPs and SNPs in the interval.
header = file_stream.readline()
for line in file_stream:
	lines = line.strip().split() 
	my_chrom = int(lines[chrom-1])
	my_pos = int(lines[pos-1])
	my_id = lines[ident-1]
	my_pval = float(lines[pval-1])
	my_zscore = float(lines[zscore-1])
	evaluate_snp(my_id,my_chrom,my_pos,my_pval, my_zscore)

#Close all filehandles.
for filename in files_dict:
	files_dict[filename].close()