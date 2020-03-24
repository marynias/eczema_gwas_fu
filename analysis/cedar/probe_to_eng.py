#! /usr/bin/env python
import os, sys, re, json
from sys import argv
#Convert illumina probe IDs in the expression file to the gene names to then ENG gene codes.

script, probes, gxp, ensembl = argv

prob_match = {}
with open (probes, 'r') as ph:
	for line in ph:
		lines = line.strip().split()
		prob_match[lines[0]] = lines[6]

with open (ensembl, 'r') as ensembl_h:
	ensembl_dict = json.load(ensembl_h)

with open (gxp, 'r') as gxph:
		header = gxph.readline()
		rest_of_file = gxph.read()
		my_ids = header.strip().split()
		new_header=""
		for my_i in my_ids:
			new_id_gene = prob_match.get(my_i, my_i)
			new_id = ensembl_dict.get(new_id_gene, new_id_gene)
			addition = new_id + " "
			new_header += addition
		print(new_header)
		print(rest_of_file)