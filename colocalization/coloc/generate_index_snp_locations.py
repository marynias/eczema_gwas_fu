#! /usr/bin/env python
#Read in a file with SNPs names and generate output file with SNP names, 
#chrosomoe and position based on Lavinia's GWAS table.

import sys

snp_list = sys.argv[1]
gwas_results = sys.argv[2]

my_snps = {}
with open (snp_list) as snp_list_h:
	for line in snp_list_h:
		lines = line.strip()
		my_snps[lines] = 1

out = open ("paternoster_2015_index_snps.txt", 'w')
with open (gwas_results) as gh:
	for line in gh:
		lines = line.strip().split()
		if lines[0] in my_snps:
			out.write("\t".join(lines[0:3]) + "\n")
out.close()