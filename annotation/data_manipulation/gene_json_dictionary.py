#! /usr/bin/env python
import os, sys, re, json
from sys import argv

script, hugo_file = argv

hugo_dict = {}
ensembl_dict = {}


def add_to_dict(l, ensembl_id):
	#Eliminate values starting with digits (they will not be IDs)
	m = re.search(r'^\d+', l.strip())
	if m is None:
		ensembl_dict[l.strip()] = ensembl_id
		hugo_dict[l.strip()] = lines[1].strip()


with open (hugo_file, 'r', encoding="utf-8") as hf:
	for line in hf:
		lines = line.strip().split("\t")
		my_eng = re.search(r"ENSG\w+", line)
		if my_eng:
			ensembl_id = my_eng.group(0)
			for l in lines:
				add_to_dict(l, ensembl_id)
			for m in lines[3].strip().split(","):
				add_to_dict(m, ensembl_id)
			for k in lines[4].strip().split(","):
				add_to_dict(k, ensembl_id)
			for n in lines[6].strip().split(","):
				add_to_dict(n, ensembl_id)

ensembl_dump = json.dumps(ensembl_dict, indent=4, sort_keys=True)
hugo_dump = json.dumps(hugo_dict, indent=4, sort_keys=True)

et = open("hugo_synonyms_ids2_filtered.ensembl.table", 'w')
ht = open("hugo_synonyms_ids2_filtered.hugo.table", 'w')
for k in sorted(ensembl_dict.keys()):
	et.write(k + "\t" + ensembl_dict[k] + "\n")
for k in sorted(hugo_dict.keys()):
	ht.write(k + "\t" + hugo_dict[k] + "\n")
et.close()
ht.close()


#eo = open("hugo_synonyms_ids2_filtered.ensembl", 'w')
#ho = open("hugo_synonyms_ids2_filtered.hugo", 'w')
#eo.write(ensembl_dump)
#ho.write(hugo_dump)
#eo.close()
#ho.close()