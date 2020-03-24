#! /usr/bin/env python
import os, sys
import argparse, gzip
from collections import defaultdict as dd

#Subset the SNPs in a given file to the given interval. 

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with SNP')
ap.add_argument('--header_tab',required=False,type=str,help='Is header present in the table? Y for yes, otherwise assumed that not')
ap.add_argument('--pos_tab',required=True,type=int,help='Column with position in the SNP table')
ap.add_argument('--int_start',required=True,type=int,help='Start of the interval to be included in the filtered output')
ap.add_argument('--int_end',required=True,type=int,help='End of the interval to be included in the filtered output')
ap.add_argument('--output',required=False,type=str,help='Name of the output file with filtered SNP results')

args = ap.parse_args()
#Initialize
tab = args.tab
header_tab = args.header_tab
pos_tab = args.pos_tab
int_start = args.int_start
int_end = args.int_end
output = args.output

out_fh = open(output, 'w')

file_stream = open(tab, 'r')
if header_tab== "Y":
	header = file_stream.readline()
	out_fh.write(header)
for line in file_stream:
	lines = line.strip().split("\t")
	if (int(lines[pos_tab-1]) >= int_start and int(lines[pos_tab-1]) <= int_end):
		out_fh.write(line)
	else:
		pass
out_fh.close()