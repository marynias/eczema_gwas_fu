#! /usr/bin/env python
from collections import defaultdict as dd
from sys import argv

script, matrixeqtl, bedfile, output = argv

reference = dd(list)
#Save chromosome, and position of each SNP in the TwinsUK dataset.
with open(bedfile) as bh:
	for line in bh:
		lines = line.strip().split()
		my_rsid = lines[3]
		my_chrom = lines[0]
		my_pos = lines[1]
		reference[my_rsid] = [my_chrom, my_pos]

oh = open(output, 'w')
with open(matrixeqtl) as mh:
		headers = mh.readline()
		for line in mh:
			lines = line.strip().split()
			my_rsid = lines[1]
			my_gene = lines[8]
			my_beta = lines[5]
			my_se = lines[6]
			my_pvalue = lines[3]
			oh.write(my_rsid + "\t" + "\t".join(reference[my_rsid]) + "\t" + my_pvalue + "\t" + my_beta + "\t" + my_se + "\t" + my_gene + "\n")

oh.close()