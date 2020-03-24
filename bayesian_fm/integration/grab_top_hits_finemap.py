#! /usr/bin/env python
import os, sys, re, glob
from collections import defaultdict as dd
from sys import argv

script, chrom, snp, interval, bf_threshold, pp_threshold  = argv

grep_seq = 'sub_finemap.sh.o' + '*'
std_output = glob.iglob(grep_seq)

#Check if regional BF is higher than the threshold.
for s in std_output:
	with open (s) as sh:
		for line in sh:
			lines=line.strip().split(":")
			if lines[0].startswith("- Log10-BF"):
				if float(lines[1]) > float(bf_threshold):
					pass
				else:
					raise ValueError ("Log10-BF for any causal SNP in the FINEMAP run for chr %s %s %s is lower than %s" % (chrom, snp, interval, bf_threshold) )

#Analyze individual SNPs and their BF
snp_filename = "chr" + chrom + "." + snp + "." + interval + "." + "snp"
with open (snp_filename) as fh:
	header = fh.readline()
	for line in fh:
		lines=line.strip().split()
		my_rsid = lines[1]
		my_pp = lines[10]
		my_logbf = lines[11]
		#Check if logBF is high enough
		if (float(my_logbf) > float(bf_threshold) or my_logbf=="inf"):
			if float(my_pp) > float(pp_threshold):
				print (my_rsid)