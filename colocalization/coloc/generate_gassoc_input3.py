#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv

#Generate gassoc input to plot the true SNP with the highest PPH4 for top gene in the plot, and for the rest of the genes plot fake SNPs in the middle of the gene with true PPh4 for a given gene.

script, top_gene, other_genes, middle, out = argv

with open (top_gene) as tg:
	line = tg.readline()
	lines = line.strip().split()
	my_top_gene = lines[4]
	my_first_line = lines[0:4]
	my_chromosome = lines[1]

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
