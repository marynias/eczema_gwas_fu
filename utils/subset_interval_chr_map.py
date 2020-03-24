#! /usr/bin/env python
import os, sys
import argparse, gzip
from collections import defaultdict as dd

#Subset the SNPs in a given file to the given interval, also considering the right chromosome, and using BED file as input 
#containing regions to be filtered. 

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with SNP')
ap.add_argument('--header_tab',required=False,type=str,help='Is header present in the table? Y for yes, otherwise assumed that not')
ap.add_argument('--bed',required=True,type=str,help='BED file with intervals to be filtered')
ap.add_argument('--pos_tab',required=True,type=int,help='Column with position in the SNP table')
ap.add_argument('--chr_tab',required=True,type=int,help='Column with chromosome in the SNP table')
ap.add_argument('--output',required=False,type=str,help='Name of the output file with filtered SNP results')

args = ap.parse_args()
#Initialize
tab = args.tab
header_tab = args.header_tab
pos_tab = args.pos_tab
chr_tab = args.chr_tab
output = args.output
bed = args.bed

to_retain = dd(list)
#Open BED file and read in positions to be retained.
with open(bed, 'r') as bed_fh:
	for line in bed_fh:
		lines = line.strip().split()
		to_retain[lines[0]].append((lines[1], lines[2], lines[3]))

out_fh = open(output, 'w')

file_stream = open(tab, 'r')
if header_tab== "Y":
	header = file_stream.readline()
	out_fh.write(header)
for line in file_stream:
	lines = line.strip().split("\t")
	if lines[chr_tab-1] in to_retain:
		for my_tuple in to_retain[lines[chr_tab-1]]:
			if (int(lines[pos_tab-1]) >= int(my_tuple[0]) and int(lines[pos_tab-1]) <= int(my_tuple[1])):
				out_fh.write("\t".join(lines) + "\t" + my_tuple[2] + "\n")
	else:
		pass
out_fh.close()