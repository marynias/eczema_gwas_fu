#! /usr/bin/env python
import os, sys
from sys import argv
from collections import defaultdict as dd
import glob


script, gene_names, snp_list, snp_name, my_interval = argv

#Read in the names of RSIDs to be extracted.
genes = set()
rsid = {}
with open (gene_names) as gn:
	for line in gn:
		lines = line.strip().split()
		genes.add(lines[0])

with open (snp_list) as sl:
	header = sl.readline()
	for line in sl:
		lines = line.strip().split()
		rsid[lines[11]] = 1

#Open relevant gene pQTL files and hunt down relevant rsids.
for gene in genes:
	grep_seq = gene + '*'
	filename = glob.iglob(grep_seq)
	out_file = snp_name + "_" + my_interval + "_" + gene + ".sun"
	out = open(out_file, "w")
	out.write("snp     effect_allele   other_allele    effect_allele_freq      beta    se      p       n    gene_name" + "\n")
	for f in filename:
		with open(f) as fh:
			header = fh.readline()
			for line in fh:
				lines = line.strip().split()
				if lines[0] in rsid:
					out.write(line.strip() + "\t" + gene + "\n")
	out.close()
