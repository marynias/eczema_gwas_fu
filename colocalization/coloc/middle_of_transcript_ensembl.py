#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv

#Generate a file with gene names and the middle of their transcript tabulated based on GTF.

script, my_gtf = argv

processed = {}
with open (my_gtf) as proces_h:
	for line in proces_h:
		if line.startswith("#"):
			pass
		else:
			lines = line.strip().split("\t")
			if lines[2] == "transcript":
				if int(lines[4]) > int(lines[3]):
					difference = float(int(lines[4]) - int(lines[3])) / 2
					middle = int(lines[3]) + difference
				elif int(lines[3]) > int(lines[4]):
					difference = float(int(lines[3]) - int(lines[4])) / 2
					middle = int(lines[4]) + difference
				annotation = lines[8].split(";")
				names = annotation[0].strip().split(" ")
				final_names = names[1].split(".")
				if final_names[0] in processed:
					pass
				else:
					print (final_names[0] + '\t' + str(int(middle)))
					processed[final_names[0]] = 1