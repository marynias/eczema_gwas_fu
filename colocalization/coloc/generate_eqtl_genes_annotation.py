#! /usr/bin/env python

import os, sys, re, argparse
from sys import argv
import gzip
from collections import defaultdict as dd

#Generate eQTL gene annotation for SNP eQTL association at a given threshold.

#output
ensemble_ids = dd(set)

script, eqtl, threshold = argv

if eqtl.endswith(".gz"):
	file_stream = gzip.open(eqtl, 'rt')
else:
	file_stream = open(eqtl, 'r')
header = file_stream.readline()
for line in file_stream:
	lines = line.strip().split("\t") 
	my_id = lines[1]
	pval = float(lines[0])
	if pval <= float(threshold):
		ensemble = lines[16].split(",")
		[ensemble_ids[my_id].add(e.strip()) for e in ensemble if e != '']
	else:
		continue

out = open(eqtl  + "_" + threshold + ".annotable", 'w')
#Print the dictionary to file.
for k in ensemble_ids:
		out.write(k + "\t" + ",".join(ensemble_ids[k]) + "\n") 

out.close()