#! /usr/bin/env python
from sys import argv
import glob

script, my_folder = argv

pops_output = glob.iglob(my_folder+"/*")


hgnc_data = dict()
#Read in Data from HGNC into dictionary
with open ("hgnc_complete_set.txt") as hgnc_set:
	header = hgnc_set.readline().strip().split("\t")
	index_s = header.index("symbol")
	index_syn = header.index("alias_symbol")
	index_l = header.index("location")
	index_e = header.index("ensembl_gene_id")
	for line in hgnc_set:
		lines=line.strip().split("\t")
		my_symbol = lines[index_s]
		my_synonyms = lines[index_syn]
		my_location = lines[index_l]
		if (len(lines) > 19 and lines[19] != "" and lines[0] != ""):  
			my_ensg = lines[index_e]
		hgnc_data[my_ensg] = [my_location, my_symbol, my_synonyms]

	#Save as dictionary with keys: ENSG ids

print("ENSG_ID\tLocation\tSymbol\tAlias_symbols\tScore")
for s in pops_output:
	with open (s) as sh:
		header = sh.readline()
		for line in sh:
			lines=line.strip().split()
			my_ensg = lines[0]
			if my_ensg in hgnc_data:
				to_print = [my_ensg, hgnc_data[my_ensg][0], hgnc_data[my_ensg][1], hgnc_data[my_ensg][2], lines[1]]
				print("\t".join(to_print))

