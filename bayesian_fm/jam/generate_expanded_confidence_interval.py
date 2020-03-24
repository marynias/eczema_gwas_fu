#! /usr/bin/env python
import os, sys, gzip
from sys import argv
from collections import defaultdict as dd

script, ld_file, pruner_file, cs_file, threshold, output = argv 

cs_loci = set()
#Read in expanded confidence interval.
with open(cs_file) as cs_fh:
	cs_l = cs_fh.readlines()
	[cs_loci.add(x.strip()) for x in cs_l]

non_selected = {}
#Dict with pos-rsid key-val pair for SNPs in the credible interval.
credible = {}
covered = {}

#Collect positions and ids of loci that have not been picked up for analysis.
#The same for ones that have been selected.
with open(pruner_file) as ph:
	header = ph.readline()
	for line in ph:
		lines = line.strip().split()
		if lines[6] == '1':
			if lines[0] in cs_loci:
				credible[lines[2]] = lines[0]
		else:
			non_selected[lines[2]] = lines[0]

open_fh = open(output, 'w')
open_fh.write("credible_set_snp" + "\t" + "pruned_snp" + "\t" + "r2" + "\n")
ld_fh = gzip.open(ld_file, 'rt')
header = ld_fh.readline()
for line in ld_fh:
	lines = line.strip().split()
	if lines[1] in credible:
		if (lines[2] in non_selected and float(lines[4]) >= float(threshold)):
			covered[credible[lines[1]]] = 1
			open_fh.write(credible[lines[1]] + "\t" + non_selected[lines[2]] + "\t" + lines[4] + "\n")
	elif lines[1] in non_selected:	
		if (lines[2] in credible and float(lines[4]) >= float(threshold)):
			covered[credible[lines[2]]] = 1
			open_fh.write(credible[lines[2]] + "\t" + non_selected[lines[1]] + "\t" + lines[4] + "\n")
			
ld_fh.close()

for k in cs_loci:
	if k in covered:
		pass
	else:
		open_fh.write(k + "\n")
		
open_fh.close()

