#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv

#Generate gassoc input to plot fake SNPs in the middle of the gene with true PPH4 for a given gene.
script, other_genes, middle, my_chromosome, out = argv

#Read in the middle position of genes and save into memory.
middle_index = {}
with open (middle) as m:
	for line in m:
		lines = line.strip().split("\t")
		middle_index[lines[0]] = lines[1]

output = open(out, 'w')
#Read in the PPH4 file and generate gassocplot input
with open (other_genes) as og:
	for line in og:
		lines = line.strip().split()
		my_middle_position = middle_index.get(lines[1], 0)
		if my_middle_position:
			output.write(lines[1] + "\t" + my_chromosome + "\t" + my_middle_position + "\t" + lines[2] + "\n")
output.close()