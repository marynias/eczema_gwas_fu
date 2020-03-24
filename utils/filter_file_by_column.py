#! /usr/bin/env python
import os, sys
import argparse, gzip
from collections import defaultdict as dd
#Filter a table by values in a given column matching a list of IDs supplied in another file.

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with IDs')
ap.add_argument('--ref',required=True,type=str,help='Input table with reference, which is going to be filtered against a list of IDs in --tab')
ap.add_argument('--header_tab',required=False,type=str,help='Is header present in the ID table? Y for yes, otherwise assumed that not')
ap.add_argument('--header_ref',required=False,type=str,help='Is header present in the ref table? Y for yes, otherwise assumed that not')
ap.add_argument('--rsid_tab',required=True,type=int,help='Column with IDs in the results table')
ap.add_argument('--rsid_ref',required=True,type=int,help='Column with IDs in the reference table')

args = ap.parse_args()
#Initialize
tab = args.tab
ref = args.ref
header_tab = args.header_tab
header_ref = args.header_ref
rsid_tab = args.rsid_tab
rsid_ref = args.rsid_ref

rsid = {}

if tab.endswith(".gz"):
	file_stream = gzip.open(tab, 'rt')
else:
	file_stream = open(tab, 'r')
if header_tab == "Y":
	header = file_stream.readline()
for line in file_stream:
	lines = line.strip().split()
	my_rsid = lines[rsid_tab-1]
	rsid[my_rsid] = 1

if ref.endswith(".gz"):
	fh2 = gzip.open(ref, 'rt')
else:
	fh2 = open(ref, 'r')

if header_ref == "Y":
	header = fh2.readline()
	#print(header.strip())
for line in fh2:
	lines = line.strip().split()
	my_rsid = lines[rsid_ref-1]
	if my_rsid in rsid:
		print(line.strip())
