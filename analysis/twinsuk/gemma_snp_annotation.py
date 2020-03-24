#! /usr/bin/env python
import os, sys
from sys import argv

#Take a given results file and harmonize the rsids between the file and reference, transfering the rsids off the given reference. 
#Check that chroms, positions and alleles match before assigning the new rsid. 
#Input accepted: BIM files. Outputs a new BIM file for the reference.

script, my_ref, my_res = argv
filename2 = os.path.splitext(my_res)

rsid = {}


#Read in rsids of interest from the bimbam file.
with open (my_res) as fh2:
	for line in fh2:
		lines = line.strip().split(",")
		rsid[lines[0]] = 1

output = filename2[0] + ".snps"
out_fh = open(output, 'w')

#Collect all the data from the reference file.
with open (my_ref) as fh1:
	for line in fh1:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		if snp_id in rsid:
			out_fh.write(snp_id + " " + pos + " " + chrom + "\n")


out_fh.close