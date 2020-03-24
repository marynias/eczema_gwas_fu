#! /usr/bin/env python
import os, sys, re
from collections import defaultdict as dd
from sys import argv

script, chrom, snp, interval, my_list, my_map, output  = argv

output_filename = os.path.basename(output)

#Obtain the position of each SNP on the given chromosome.
positions = dict()
with open (my_map) as map_fh:
	for line in map_fh:
		lines=line.strip().split()
		if lines[1] == chrom:
			positions[lines[0]] = lines[2]
selected = dict()
#Read in the list of selected SNPs.
with open (my_list) as mh:
	for line in mh:
		selected[line.strip()] = positions[line.strip()]

out = open (output, 'w')
#Sort by selected SNPs position
my_temp = sorted(selected.items(), key=lambda x:x[1])

out.write("run")
for (k, v) in my_temp:
	out.write("\t" + k)
out.write("\n")

seen = dict()

#Analyze individual SNPs and their BF
snp_filename = "chr" + chrom + "." + snp + "." + interval + "_full_table.txt" 
with open (snp_filename) as fh:
	header = fh.readline()
	for line in fh:
		lines=line.strip().split()
		my_rsid = lines[0]
		my_pp = lines[1]
		if my_rsid in selected:
			seen[my_rsid] = my_pp

out.write(output_filename)
for (k, v) in my_temp:
	if k in seen:
		out.write("\t" + seen[k])
	else:
		out.write("\t" + "NA")
out.write("\n")
out.close()