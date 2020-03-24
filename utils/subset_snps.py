#! /usr/bin/env python
import os, sys
import argparse, gzip
from collections import defaultdict as dd

#Subset the SNPs in a given file to the given list of SNPs.

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with SNP')
ap.add_argument('--snp_list',required=True,type=str,help='File with list of SNPs to be filtered')
ap.add_argument('--header_tab',required=False,type=str,help='Is header present in the table? Y for yes, otherwise assumed that not')
ap.add_argument('--pos_tab',required=True,type=int,help='Column with SNP ID in the SNP table')
ap.add_argument('--output',required=False,type=str,help='Name of the output file with filtered SNP results')

args = ap.parse_args()
#Initialize
tab = args.tab
header_tab = args.header_tab
snp_list = args.snp_list
pos_tab = args.pos_tab
output = args.output

to_retain = {}
#Open BED file and read in positions to be retained.
with open(snp_list, 'r') as snp_list_fh:
	for line in snp_list_fh:
		lines = line.strip().split()
		to_retain[lines[0]] = 1

out_fh = open(output, 'w')

file_stream = open(tab, 'r')
if header_tab== "Y":
	header = file_stream.readline()
	out_fh.write(header)
for line in file_stream:
	lines = line.strip().split("\t")
	if lines[pos_tab-1] in to_retain:
		out_fh.write(line)
	else:
		pass
out_fh.close()