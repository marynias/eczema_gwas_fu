#! /usr/bin/env python
import os, sys, re, gzip
from sys import argv
from collections import defaultdict as dd
#Extract Blueprint data for relevant genes and tissues.
script, my_blueprint, my_list, gn, tissue = argv

#Read in the list of genes and lead SNPs to look for.
genes = set()
my_rsid_match = dd(set)
with open(my_list) as my_list_fh:
	for line in my_list_fh:
		lines = line.strip().split()
		genes.add(lines[0])
		my_rsid_match[lines[1]].add(lines[0])
#Read in gene names.
gene_names = dd()
with open(gn) as gene_names_fh:
	for line in gene_names_fh:
		lines = line.strip().split()
		gene_names[lines[0]] = lines[2]

#Iterate over Blueprint file.
saved = dd(set)
if my_blueprint.endswith(".gz"):
	file_stream = gzip.open(my_blueprint, 'rt')
else:
	file_stream = open(my_blueprint, 'r')

header = "\t".join(["chr:pos_ref_alt", "rsid", "phenotypeID", "p.value", "beta", "Bonferroni.p.value", "FDR", "alt_allele_frequency", "std.error_of_beta"   ])

for line in file_stream:
	lines = line.strip().split() 
	(my_gene, my_suffix) = lines[2].strip().split(".")
	if my_gene in genes:
		saved[my_gene].add("\t".join(lines))
file_stream.close()

#Print matching lines:
for rsid in my_rsid_match:
	for my_locus in my_rsid_match[rsid]:
		if my_locus in saved:
			gene_name = gene_names.get(my_locus, my_locus)
			output = rsid + "_" + gene_name + "_3Mbp_" + tissue + ".blueprint"
			temp_fh = open(output, 'w')
			temp_fh.write(header + "\t" + "external_gene_name" + "\n")
			for my_lines in saved[my_locus]:
				temp_fh.write(my_lines.strip() + "\t" + gene_name + "\n")
			temp_fh.close()

