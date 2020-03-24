#! /usr/bin/env python

import os, sys, re
from sys import argv
import gzip
from collections import defaultdict as dd

script, significant, eqtl = argv

my_targets = dd()
file_stream = open(eqtl, 'r')
for line in file_stream:
	lines = line.strip().split()  
	my_targets[lines[0]] = lines[1]
file_stream.close()

#Read in a list of significant SNPs, annotate and print out the results.
#First check if file empty - if so, terminate the script early.
if os.stat(significant).st_size == 0:
	sys.exit()

out = open(significant + ".annotated", 'w')
out.write("snp" + "\t" + "pvalues.gwas" + "\t" + "MAF.gwas" + "\t" + "pvalues.eqtl" + "\t" + "SNP.PP.H4" + "\t" + "Ensemble_ID" + "\n")

with open (significant) as proces_h:
	for line in proces_h:
		lines = line.strip().split()  
		#Values: P-value in df1 (GWAS), MAF in df1 (GWAS), pvalue in df2 (eQTL), SNP.PP.H4
		#Retrieve relevant gene annotation IDs.
		if lines[0] in my_targets:
			targets = my_targets[lines[0]].split(",")
			[out.write(lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\t" + lines[7] + "\t" \
			+ lines[14] + "\t" + t + "\n") for t in targets]
		else:
			pass

out.close()