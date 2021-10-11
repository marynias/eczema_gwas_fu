#! /usr/bin/env python
import os, sys
from sys import argv
from collections import defaultdict as dd
from Bio import SeqIO

script, template, vcf, reference = argv

#Load in the variant file
variant_vcf = dd(lambda: dd())
with open (vcf) as vcf_h:
	for line in vcf_h:
		if line.startswith("#"):
			pass
		else:
			lines = line.strip().split()
			chrom = lines[0][3:]
			pos = lines[1]
 			ref = lines[3]
 			alt = lines[4]
 			genotypes = lines[9].split(":")
 			genotype = genotypes[0]
 			if genotype == "0/1":
 				saved_genotype = ref + alt
 			elif genotype == "1/0":
 				saved_genotype = ref + alt
 			elif genotype == "1/1":
 				saved_genotype = alt + alt
 			else:
 				saved_genotype = ref + ref
 			variant_vcf[chrom][pos] = saved_genotype

#Load in the reference file
my_reference = dd()
records = list(SeqIO.parse(reference, "fasta"))
for r in records:
	chrom = r.id
	my_chrom = chrom[3:]
	my_reference[my_chrom] = r.seq.upper()

#Load in the template file.
var = 0
with open (template) as fh1:
	for line in fh1:
		if line.startswith("#"):
			lines = line.strip().split()
			print ("\t".join(lines))
		else:
			lines = line.strip().split()
			if (lines[1] in variant_vcf and lines[2] in variant_vcf[lines[1]]):
				print (lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\t" + variant_vcf[lines[1]][lines[2]])
				var += 1
				#print ("var")
			elif lines[1] in my_reference:
				 print (lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\t" + my_reference[lines[1]][int(lines[2])-1] + my_reference[lines[1]][int(lines[2])-1])
				 #print("genome")
			else:
				print ("\t".join(lines))
	#print (var)
