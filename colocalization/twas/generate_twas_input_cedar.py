#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict as dd

script, ensembl, rpkm, fam, output_prefix = argv

#Read in Ensembl gene names.
ensembl_ids = {}
with open (ensembl) as ensembl_h:
	for line in ensembl_h:
		lines = line.strip().split()
		ensembl_ids[lines[1]] = lines[0]

#Read in RPKM values and create appropriate output files.
rpkm_values = {}
files_dict = {}
with open (rpkm) as rpkm_values_h:
	header = rpkm_values_h.readline()
	headers = header.strip().split()
	for line in rpkm_values_h:
		lines = line.strip().split()
		gene_id = lines[0]
		gene_name = ensembl_ids.get(gene_id, gene_id)
		rpkm_values[gene_name] = lines[0:]
		filename = output_prefix + "_" + gene_name + ".fam"
		files_dict.setdefault(filename, open(filename,'w')) 

#Read in the fam file.
with open (fam) as fam_h:
	for line in fam_h:
		lines = line.strip().split()
		my_entry = " ".join(lines[0:5])
		my_ind = lines[0]
		for gene_name in rpkm_values:
			filename = output_prefix + "_" + gene_name + ".fam"
			files_dict[filename].write(my_entry)
			#Write the correct RPKM value.
			my_index = headers.index(my_ind)
			print (fam)
			print (gene_name)
			print (my_ind)
			print (my_index)
			my_rpkm = rpkm_values[gene_name][my_index]
			files_dict[filename].write(" " + my_rpkm + "\n")

#Close all filehandles.
for filename in files_dict:
	files_dict[filename].close()