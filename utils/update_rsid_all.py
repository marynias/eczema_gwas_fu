#! /usr/bin/env python
import os, sys
import argparse
from collections import defaultdict as dd

#Update a given results table with rsids from a reference contained in a bim file, 
#In the case, when the input SNP not found in reference, print out its original name.

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with results')
ap.add_argument('--bim',required=True,type=str,help='Input BIM with reference rsids')
ap.add_argument('--head',required=False,type=str,help='Is header present? Y for yes, otherwise assumed that not')
ap.add_argument('--chrom',required=True,type=int,help='Column with chromosome number in the input table')
ap.add_argument('--pos',required=True,type=int,help='Column with position number in the input table')
ap.add_argument('--ref',required=True,type=int,help='Column with reference allele in the input table')
ap.add_argument('--alt',required=True,type=int,help='Column with alternative allele in the input table')
ap.add_argument('--rsid',required=False,type=int,default=-999,help='Column with rsid id in the input table, if present')

args = ap.parse_args()
#Initialize
tab = args.tab
bim = args.bim
chrom = args.chrom
pos = args.pos
ref = args.ref
alt = args.alt
my_rsid = args.rsid
head = args.head

loci = dd(lambda: dd(set))
rsid = {}

#Collect all the data from the reference file.
with open (bim) as fh1:
	for line in fh1:
		lines = line.strip().split()
		my_chrom = lines[0]
		my_pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		loci[my_chrom][my_pos].add(snp_id)
		rsid[snp_id] = [allele_1, allele_2]

with open (tab) as fh2:
	if head == "Y":
		if my_rsid == -999:
			header = fh2.readline()
			print (header.strip() + "\t" + "RSID")
		else:
			header = fh2.readline()
			print (header.strip())
			
	for line in fh2:
		lines = line.strip().split()
		if my_rsid == -999:
			lines.append("")
			my_r = -1	
		else:
			my_r = my_rsid - 1
		my_chrom = lines[chrom-1]
		my_pos = lines[pos-1]
		allele1_1 = lines[ref-1]
		allele1_2 = lines[alt-1]
		if my_chrom in loci:
			if my_pos in loci[my_chrom]:
				my_rsids = loci[my_chrom][my_pos]
				for r in my_rsids:
					allele2_1 = rsid[r][0]
					allele2_2 = rsid[r][1]
					if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
						lines[my_r] = r
						print ("\t".join(lines))
					elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
						lines[my_r] = r
						print ("\t".join(lines))
					else:
						print ("\t".join(lines))
			else:
				print ("\t".join(lines))

		else:
			print ("\t".join(lines))
