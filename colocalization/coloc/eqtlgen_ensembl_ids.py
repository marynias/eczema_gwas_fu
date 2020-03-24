#! /usr/bin/env python

import os, sys, re
import gzip
from collections import defaultdict as dd
from sys import argv

#Substitute transcript IDS in eqtlgen files with gene names within 3MBbp of index SNP in eczema GWAS.
#Print a subset of eQTLgen dataset containing only substituted lines.

script,ensembl,eqtlgen,output = argv

ensembl_ids = {}
proc = dd(lambda: dd(int))

with open (ensembl) as ensembl_h:
	for line in ensembl_h:
		lines = line.strip().split()  
		ensembl_ids[lines[1]] = lines[0]

#Dictionary holding of SNP-gene pairs processed.

out = open(output, 'w')
if eqtlgen.endswith(".gz"):
	eqtlgen_h = gzip.open(eqtlgen, 'rt')
else:
	eqtlgen_h = open(eqtlgen, 'r')
header = eqtlgen_h.readline()
out.write(header)
for line in eqtlgen_h:
	lines = line.strip().split("\t")
	#print(lines)
	ensembl_transcript = lines[16]
	my_snp = lines[1]
	ensembl_transcripts = lines[16].strip().split(",")
	ensembl_transcripts = [x.strip() for x in ensembl_transcripts]
	for e in ensembl_transcripts:
		#Check that transcript is of interest.
		if e in ensembl_ids:
			my_gene = ensembl_ids[e]
			if my_snp in proc:
				if my_gene in proc[my_snp]:
					pass
				else:
					lines[16] = my_gene
					out.write("\t".join(lines) + "\n")
					proc[my_snp][my_gene] = 1
			else:
				lines[16] = my_gene
				out.write("\t".join(lines) + "\n")
				proc[my_snp][my_gene] = 1

eqtlgen_h.close()
out.close()