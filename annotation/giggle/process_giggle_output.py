#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv
from collections import defaultdict as dd

#Modify giggle output to print the query rsid, chrom and posiition as the first 3 entry on each line followed by the hit. 

script, giggle, output = argv

if giggle.endswith(".gz"):
	giggle_fh = gzip.open(giggle, 'rt')
else:
	giggle_fh = open(giggle, 'r')

out = open(output, 'w')

for line in giggle_fh:
	if line.startswith("##"):
		lines = line.strip().split()	
		my_rsid = lines[3]
		my_chromosome = lines[0].strip("##")
		my_position_start = lines[1] 
		my_position_end = lines[2] 
	else:
		out.write(my_chromosome + "\t" + my_position_start + "\t"
		+ my_position_end + "\t" + my_rsid + "\t" + line)

giggle_fh.close()
out.close()
