#! /usr/bin/env python
import os, sys
from sys import argv
from collections import defaultdict as dd

#Take a given results file and harmonize the rsids between the file and reference, transfering the rsids off the given reference. 
#Check that chroms, positions and alleles match before assigning the new rsid. 
#Input accepted: BIM files. Outputs a new BIM file for the reference.

script, my_ref, my_res = argv
filename1 = os.path.basename(my_ref)
filename2 = os.path.basename(my_res)

loci = dd(lambda: dd(set))
rsid = {}

#Collect all the data from the reference file.
with open (my_ref) as fh1:
	for line in fh1:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		loci[chrom][pos].add(snp_id)
		rsid[snp_id] = [allele_1, allele_2]

with open (my_res) as fh2:
	for line in fh2:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		allele1_1 = lines[4]
		allele1_2 = lines[5]
		if chrom in loci:
			if pos in loci[chrom]:
				my_rsids = loci[chrom][pos]
				for r in my_rsids:
					allele2_1 = rsid[r][0]
					allele2_2 = rsid[r][1]
					if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
						print (chrom + "\t" + r + "\t" + "0" + "\t" + pos + "\t" + allele1_1 + "\t" + allele1_2)
					elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
						print (chrom + "\t" + r + "\t" + "0" + "\t" + pos + "\t" + allele1_1 + "\t" + allele1_2)
					else:
						continue
