#! /usr/bin/env python
import os, sys
import argparse, gzip
from collections import defaultdict as dd

#Harmonize the sign of beta in two files: make sure that beta is given relative to the same reference allele.

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with results')
ap.add_argument('--ref',required=True,type=str,help='Input table with reference allele, relative to which beta in --tab is going to adjusted to')
ap.add_argument('--header_tab',required=False,type=str,help='Is header present in the results table? Y for yes, otherwise assumed that not')
ap.add_argument('--header_ref',required=False,type=str,help='Is header present in the ref table? Y for yes, otherwise assumed that not')
ap.add_argument('--rsid_tab',required=True,type=int,help='Column with rsid in the results table')
ap.add_argument('--rsid_ref',required=True,type=int,help='Column with rsid in the reference table')
ap.add_argument('--effect_tab',required=True,type=int,help='Column with effect allele in the results table')
ap.add_argument('--alt_tab',required=True,type=int,help='Column with the alternative allele in the results table')
ap.add_argument('--effect_ref',required=True,type=int,help='Column with effect allele in the reference table')
ap.add_argument('--beta_tab',required=True,type=int,help='Column with beta value in the results table')
ap.add_argument('--zscore_tab',required=False,type=int,help='Column with zscore values in the results table, if present')
ap.add_argument('--out',required=True,type=str,help='Name of output file')

args = ap.parse_args()
#Initialize
tab = args.tab
ref = args.ref
header_tab = args.header_tab
header_ref = args.header_ref
rsid_tab = args.rsid_tab
rsid_ref = args.rsid_ref
effect_tab = args.effect_tab
alt_tab = args.alt_tab
effect_ref = args.effect_ref
beta_tab = args.beta_tab
zscore_tab = args.zscore_tab
out = args.out

rsid = {}

if ref.endswith(".gz"):
	file_stream = gzip.open(ref, 'rt')
else:
	file_stream = open(ref, 'r')

if header_ref == "Y":
	header = file_stream.readline()
for line in file_stream:
	lines = line.strip().split("\t")
	my_rsid = lines[rsid_ref-1]
	effect_allele = lines[effect_ref-1]
	if my_rsid in rsid:
		if effect_allele == rsid[my_rsid]:
			continue
		else:
			print ("Current effect allele: %s" % effect_allele)
			print ("Previous effect allele: %s" % rsid[my_rsid])
			print ("Current snp: %s" % my_rsid)
			raise ValueError ("Multiple effect alleles for the same SNP in the reference.")
	else:
		rsid[my_rsid] = effect_allele

out_fh = open(out, 'w')
with open (tab) as fh2:
	if header_tab == "Y":
		header = fh2.readline()
		out_fh.write(header)
	for line in fh2:
		lines = line.strip().split()
		my_rsid = lines[rsid_tab-1]
		effect_allele = lines[effect_tab-1]
		alt_allele = lines[alt_tab-1]
		if my_rsid in rsid:
			if rsid[my_rsid] == effect_allele:
				out_fh.write("\t".join(lines) + "\n")
			else:
				current_beta = lines[beta_tab-1]
				new_beta = -float(current_beta)
				lines[beta_tab-1] = new_beta
				lines[effect_tab-1] = alt_allele
				lines[alt_tab-1] = effect_allele
				if zscore_tab:
					current_z = lines[zscore_tab-1]
					new_z = -float(current_z)
					lines[zscore_tab-1] = new_z
				print ("Allele effect direction changed for snp  %s" % my_rsid)
				lines = [str(x) for x in lines]
				out_fh.write("\t".join(lines) + "\n")
		else:
			print("Variant %s not present in the reference file" %my_rsid)

out_fh.close()