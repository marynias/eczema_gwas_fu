#! /usr/bin/env python
import os, sys, re
from collections import defaultdict as dd
from sys import argv

script, chrom, snp, interval, bf_threshold, pp_threshold = argv

#Check if regional BF is higher than the threshold.
s = "Log.BayesFactor"

with open (s) as sh:
	all_lines = sh.readlines()
	bf_line = all_lines[0].strip()
	if (float(bf_line) >float(bf_threshold) or bf_line == "Inf"):
		pass
	else:
		raise ValueError ("Log10-BF for any causal SNP in the Paintor run for chr %s %s %s is lower than %s" % (chrom, snp, interval, bf_threshold) )



snp_filename = "chr" + chrom + "." + snp + "." + interval + ".results" 
with open (snp_filename) as fh:
	header = fh.readline()
	for line in fh:
		lines=line.strip().split()
		my_rsid = lines[2]
		my_pp = lines[6]
		if float(my_pp) > float(pp_threshold):
			print (my_rsid)