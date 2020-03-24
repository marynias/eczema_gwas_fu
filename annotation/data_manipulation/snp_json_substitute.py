#! /usr/bin/env python
import os, sys, json, re
import argparse
from collections import defaultdict as dd

ap = argparse.ArgumentParser()

ap.add_argument('--tab',required=True,type=str,help='Input table with results')
ap.add_argument('--db',required=True,type=str,help='Input JSON dictionary matching aliases to dbSNP reference database / HUGO gene ref database')
ap.add_argument('--head',required=False,type=str,help='Is header present? Y for yes, otherwise assumed that not')
ap.add_argument('--rsid',required=True,type=int,default=-999,help='Column number with rsid id / gene_id in the input table')
ap.add_argument('--delim',required=True,type=str,default=-999,help='Delimiter string used in the table. A - \t, B - space C-comma')

args = ap.parse_args()
#Initialize
tab = args.tab
db = args.db
my_rsid = args.rsid
head = args.head
my_delim = args.delim
#Load in dbSNP mappings.
with open (db, 'r') as db_h:
	db_dict = json.load(db_h)
with open (tab) as fh2:
	if head == "Y":
			header = fh2.readline()
			print (header.strip())
	for line in fh2:
		if my_delim == 'A':
			lines = line.strip().split("\t")
		elif my_delim == 'B':
			lines = line.strip().split()
		elif my_delim == 'C':
			lines = line.strip().split(",")
		my_replacement = db_dict.get(lines[my_rsid-1], lines[my_rsid-1])
		lines[my_rsid-1] = my_replacement
		print ("\t".join(lines))
