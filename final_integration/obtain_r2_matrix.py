#! /usr/bin/env python
import os, sys
import gzip
from sys import argv
from collections import defaultdict as dd
from math import sqrt
#Obtain R matrix based on 1k EUR results. Compare SNPs by position.

script, rsid_list, ld_file, out = argv
#Read in positions of interest.
rsids = dict()
rsid_set = set()
with open (rsid_list) as fh1:
	header = fh1.readline()
	for line in fh1:
		lines = line.strip().split()
		rsids[lines[2]] = lines[0]
		rsid_set.add(lines[2])

sorted_by_pos = sorted(rsids.keys(), key=int)
rsids_by_pos = [rsids[k] for k in sorted_by_pos]

if ld_file.endswith(".gz"):
	file_stream = gzip.open(ld_file, 'rt')
else:
	file_stream = open(ld_file, 'r')
saved_correlation = dd(lambda: dd(lambda: 'NA'))
header = file_stream.readline()
for line in file_stream:
	lines = line.strip().split()
	if lines[1] in rsid_set:
		if lines[2] in rsid_set:
			r2 = float(lines[4])
			r = sqrt(r2)
			saved_correlation[lines[1]][lines[2]] = str(r)
			saved_correlation[lines[2]][lines[1]] = str(r)
file_stream.close()

out_fh = open(out, 'w')
out_fh.write("id" + "\t" + "\t".join(rsids_by_pos) + "\n")
for my_value in sorted_by_pos:
	current_line = ""
	current_line += rsids[my_value]
	current_line += "\t"
	for my_value2 in sorted_by_pos:
		if my_value == my_value2:
			r = "1"
			current_line += r
			current_line += "\t"
		else:
			r = saved_correlation[my_value][my_value2]
			current_line += r
			current_line += "\t"
	out_fh.write(current_line + "\n")
out_fh.close()