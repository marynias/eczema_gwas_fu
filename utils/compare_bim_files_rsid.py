#! /usr/bin/env python

#Compare two files by rsids.
#Allele values have to be the same case (upper or lower) in both files being compared.

import os, sys
from sys import argv

script, file1, file2 = argv

loci = {}
with open (file1) as fh1:
	for line in fh1:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		rsid[snp_id] = [chrom, pos, allele_1, allele_2]

with open (file2) as fh2:
	for line in fh2:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		if snp_id in rsid:
			if chrom == rsid[snp_id][0]:
				if pos == rsid[snp_id][1]:
					if (allele_1 == rsid[snp_id][2] and allele_2 == rsid[snp_id][3]):
						continue
					elif (allele_1 == rsid[snp_id][3] and allele_2 == rsid[snp_id][2]):
						#Potentially develop this to output flipped alleles
						continue
					else:
						print ("Warning! Variant %s not matched by alleles" % snp_id)

				else:
					print ("Warning! Variant %s not matched by position" % snp_id)
			else:
				print ("Warning! Variant %s not matched by chromosome" % snp_id) 
