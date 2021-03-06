#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv

#Generate a file with gene Ensembl ID  and the start and end of their transcript tabulated based on GTF.

script, my_gtf = argv

processed = {}
with open (my_gtf) as proces_h:
	for line in proces_h:
		if line.startswith("#"):
			pass
		else:
			lines = line.strip().split("\t")
			if lines[2] == "gene":
				my_chrom = lines[0]
				if int(lines[4]) > int(lines[3]):
					start = int(lines[3]) 
					end = int(lines[4]) 
				elif int(lines[3]) > int(lines[4]):
					end = int(lines[3]) 
					start = int(lines[4]) 
				annotation = lines[8].split(";")
				names = annotation[0].strip().split(" ")
				ensembl_id = names[1].split(".")
				if ensembl_id[0] in processed:
					pass
				else:
					print (ensembl_id[0] + '\t' + str(my_chrom) + '\t' + str(start) + '\t' + str(end))
					processed[ensembl_id[0]] = 1