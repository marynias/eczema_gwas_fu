#! /usr/bin/env python
import os, sys, re, json
from sys import argv
from collections import defaultdict as dd

script, eqtls, hugo, suffix = argv

#Read in Ensembl to HUGO mappings.
with open (hugo, 'r') as hugo_h:
	hugo_dict = json.load(hugo_h)

my_lines =  dd(list)
with open (eqtls, 'r') as eqtls_h:
	for line in eqtls_h:
		lines = line.strip().split()
		my_replacement = hugo_dict.get(lines[0], lines[0])
		lines[0] = my_replacement
		line = "\t".join(lines)
		my_lines[lines[0]].append(line)

#Print output per each gene.
for gene in my_lines:
	output = gene + "_" + suffix
	current_fh = open (output, 'w')
	for l in my_lines[gene]:
		current_fh.write(l + "\n")
	current_fh.close()