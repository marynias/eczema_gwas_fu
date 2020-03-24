#! /usr/bin/env python
import os, sys, re
from collections import defaultdict as dd
from sys import argv

script, chrom, snp, interval, bf_threshold, pp_threshold = argv

#Check if regional BF is higher than the threshold.
s = "chr" + chrom + "." + snp + "." + interval + "_region_bf.txt" 
with open (s) as sh:
	all_lines = sh.readlines()
	bf_line = all_lines[3].split()
	if (float(bf_line[1]) > bf_threshold or bf_line[1] == "Inf"):
		pass
	else:
		raise ValueError ("Log10-BF for any causal SNP in the JAM run for chr %s %s %s is lower than %s" % (chrom, snp, interval, bf_threshold) )

selected = {}

snp_filename = "chr" + chrom + "." + snp + "." + interval + "_full_table.txt" 
with open (snp_filename) as fh:
	header = fh.readline()
	for line in fh:
		lines=line.strip().split()
		my_rsid = lines[0]
		my_pp = lines[1]
		my_bf = lines[8]
		if my_bf != "NA":
			if (float(my_bf) > float(bf_threshold) or my_bf=="Inf"):
				if float(my_pp) > float(pp_threshold):
					print (my_rsid)
					selected[my_rsid] = 1

#Check extended credible interval for any SNPs which are in high LD with selected SNPs above
extended_filename = "chr" + chrom + "." + snp + "." + interval + "_credible_set_extended.txt" 
with open (extended_filename) as eh:
	header = eh.readline()
	for line in eh:
		lines=line.strip().split()
		if len(lines) > 1:
			if lines[0] in selected:
				pass
				#Add additional rsid in high LD with our target from selected
				#print(lines[1])
