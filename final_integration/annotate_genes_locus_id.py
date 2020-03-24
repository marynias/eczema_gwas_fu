#! /usr/bin/env python

import os, sys
from collections import defaultdict as dd
import glob

lookups = glob.glob('/panfs/panasas01/sscm/qh18484/final_integration/lookups/processed/*processed')

lookup_pairs = dd(set)

#Read in all the gene-locus ID pairs found in among the lookups.
for l in lookups:
	filehandle = open(l, 'r')
	header = filehandle.readline()
	for line in filehandle:
		lines = line.strip().split("\t")
		if len(lines) > 5:
			lookup_pairs[lines[4].strip()].add(lines[2].strip())
	filehandle.close()

#for key in lookup_pairs:
#	print (key + "\t" + ",".join(list(lookup_pairs[key])))

genes = glob.glob('/panfs/panasas01/sscm/qh18484/final_integration/gene/processed/*processed')

#Read in all the gene-locus ID pairs found in among the lookups.
for g in genes:
	my_filename = os.path.basename(g)
	#print(my_filename)
	output = os.path.splitext(my_filename)[0] + ".index_locus.processed"
	filehandle = open(g, 'r')
	my_out = open(output, 'w')
	header = filehandle.readline()
	my_out.write(header)
	for line in filehandle:
		lines = line.strip().split("\t")
		my_gene = lines[4]
		#Print out the same line but with added index RSIDs.
		if my_gene in lookup_pairs:
			if lines[2] == 'NA':
				my_rsids = lookup_pairs[my_gene]
				for r in my_rsids:
					lines[2] = r
					my_out.write("\t".join(lines) + "\n")
			else:
				my_out.write("\t".join(lines) + "\n")
		else:
			my_out.write("\t".join(lines) + "\n")
	filehandle.close()
	my_out.close()