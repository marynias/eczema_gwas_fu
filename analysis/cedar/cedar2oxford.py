#! /usr/bin/env python
import os, sys
from sys import argv
from collections import defaultdict as dd

#Convert Dosage MatrixEQTL format to Oxford Geno.
#Less than 1 - consider for heterozygotes and null homozygote
#More than 1 - consider for heterozygote and double homozygote

script, cedar_markers, cedar_dosage = argv

snps = {}
with open (cedar_markers) as cedar_markers_h:
	header = cedar_markers_h.readline()
	for line in cedar_markers_h:
		lines = line.strip().split()
		details = lines[0].split("_")
		details.append(lines[1])
		#Following elements: chromosome, pos, ref_allele, alt_allele, rsid

		snps[lines[0]] = details
files_dict = {}
for a in range (1,23):
	filename =  "cedar.chr" + str(a) + ".gen"
	files_dict.setdefault(filename, open(filename,'w')) 

with open (cedar_dosage) as cedar_dosage_h:
	header = cedar_dosage_h.readline()
	for line in cedar_dosage_h:
		lines = line.strip().split()
		details = lines[0].split("_")
		ref_allele = details[2]
		alt_allele = details[3]
		no_of_sample = len(lines)-1
		lines[1:]=[float(x) for x in lines[1:]]
		sum_dosage = sum(lines[1:])
		if sum_dosage > no_of_sample:
#Means we are counting the dosage of major allele, but we need to reverse it to use minor allele as reference
#Have to make sure to output the probability of genotype, with minor allele being the first listed, i.e. allele A.
			(ref_allele, alt_allele) = (alt_allele, ref_allele)
			lines[1:]=[2 - x for x in lines[1:]]
		#Print info about locus
		my_info = snps[lines[0]]
		chrom = my_info[0]
		rsid = my_info[4]
		pos = my_info[1]
		filename =  "cedar.chr" + str(chrom) + ".gen"
		files_dict[filename].write(str(chrom) + " " + rsid + " " + str(pos) + " " + ref_allele + " " + alt_allele + " ")
		#For each dosage value, produce probabilities for 3 genotypes:
		for values in lines[1:]:
			if values == 0:
				files_dict[filename].write("0" + " " + "0" + " " + "1" + " ")
			elif values < 1:
				files_dict[filename].write("0" + " " + str(values) + " " + str(1-values) + " ")
			elif values == 1:
				files_dict[filename].write("0" + " " + "1" + " " + "0" + " ")
			elif values < 2:	
				files_dict[filename].write(str(values - 1) + " " + str(2 - values) + " " + "0" + " ")
			elif values == 2:	
				files_dict[filename].write("1" + " " + "0" + " " + "0" + " ")
		files_dict[filename].write("\n")	
#Close all filehandles.
for filename in files_dict:
	files_dict[filename].close()