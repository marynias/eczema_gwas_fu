#! /usr/bin/env python
import os, sys, re, json
from sys import argv

script, snp_file = argv

rsid_dict = {}

#Create a dictionary with keys = dbsnp current ref rsid, value = rsid alias

with open (snp_file, 'r') as sf:
	header = sf.readline()
	for line in sf:
		lines = line.strip().split()
		rsid_dict[lines[1]] = lines[0]

with open('rsid_synonyms.txt', 'w') as f:
        json.dump(rsid_dict, f, indent=4, sort_keys=True)