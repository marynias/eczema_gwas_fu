#! /usr/bin/env python
import os, sys, json, re
import argparse, gzip
from collections import defaultdict as dd

ap = argparse.ArgumentParser()

ap.add_argument('--tab',required=True,type=str,help='Input table with results.')
ap.add_argument('--queries',required=True,type=str,help='Input BED-like file with input SNPs and index SNPs to whose locus they belong')
ap.add_argument('--head',required=False,type=int,default=0,help='Line at which header starts; all previous lines ignored')
ap.add_argument('--my_id',required=True,type=int,help='Column number with rsid id / gene_id in the input table')
ap.add_argument('--delim',required=True,type=str,help='Delimiter string used in the table.')
ap.add_argument('--delim_c',required=True,type=str,help='Delimiter string used within cells.')

args = ap.parse_args()
#Initialize
my_table = args.tab
queries = args.queries
my_rsid = args.my_id
head = args.head
my_delim = args.delim
my_delim_c = args.delim_c

queries_index = dd(set)


#Key - SNP id, value = set of index SNPs matching the SNP id.
with open (queries, 'r') as queries_h:
	for line in queries_h:
		lines = line.strip().split()
		queries_index[lines[3]].add(lines[4])
with open (my_table, 'r') as my_table_h:
	if head:
		for line in range(head-1):
			my_table_h.readline()
		header = my_table_h.readline().strip().split()
		header.append("Index_SNP")
		print ("\t".join(header))		
	for line in my_table_h:
		lines = re.split(my_delim,line.strip())
		lines.append("Index_SNP")
		values = re.split(my_delim_c,lines[my_rsid-1])
		for v in values:
			if v.strip() in queries_index:
				matching_index_snps = queries_index[v.strip()]
				#Print one line for one matching SNP and its matching index SNP
				for m in matching_index_snps:
					lines[my_rsid-1] = v.strip()
					lines[-1] = m
					print("\t".join(lines))
