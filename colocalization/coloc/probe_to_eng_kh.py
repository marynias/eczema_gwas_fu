#! /usr/bin/env python
import os, sys, re, json
from sys import argv
#Convert illumina probe IDs in the expression file to the gene names to then ENG gene codes.

script, probes, gxp, ensembl = argv

prob_match = {}
with open (probes, 'r') as ph:
	for line in ph:
		lines = line.strip().split()
		prob_match[lines[0]] = lines[3]

with open (ensembl, 'r') as ensembl_h:
	ensembl_dict = json.load(ensembl_h)

with open (gxp, 'r') as gxph:
		for line in gxph:
			lines = line.strip().split()
			#Get Gene name.
			gene_name = prob_match.get(lines[0], lines[0])
			my_replacement = ensembl_dict.get(gene_name, gene_name)
			lines[0] = my_replacement
			print ("\t".join(lines))