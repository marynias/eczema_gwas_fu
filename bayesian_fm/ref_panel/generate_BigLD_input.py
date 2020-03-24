#! /usr/bin/env python
import os, sys, re
from collections import defaultdict as dd
from sys import argv
import gzip

script, my_vcf, out = argv

out_matrix = open(out + ".geno", 'w')
out_beta = open(out + ".info", 'w')
out_beta.write("rsID" + "\t" + "bp" + "\n")

processed = {}

def print_bigld_output(my_rsid, rows):
	out_beta.write(my_rsid + "\t" + rows[rsid_pos] + "\n")
	out_matrix.write(my_rsid)
	count_val = "1"
	for indin in rows[9:]:
		genotypes = indin.split("|")
		my_count = genotypes.count(count_val)
		out_matrix.write("\t" + str(my_count))
	out_matrix.write("\n")

#Iterate through VCF files to print genotype matrix and beta vector.
if my_vcf.endswith(".gz"):
	file_stream = gzip.open(my_vcf, 'rt')
else:
	file_stream = open(my_vcf, 'r')
for lines in file_stream:
    if(lines[0:2] != "##"):
        if (lines[0:6] == "#CHROM"):
            header = lines.strip().split()
            A0_index = header.index("REF")
            A1_index = header.index("ALT")
            rsid_index = header.index("ID")
            rsid_pos = header.index("POS")
            out_matrix.write(" " + "\t" + "\t".join(header[9:]) + "\n")
        else:
            rows = lines.strip().split()
            my_rsid = rows[rsid_index]
            vcf_A0 = rows[A0_index]
            vcf_A1 = rows[A1_index]
            if len(vcf_A0) == len(vcf_A1):
                if my_rsid in processed:
                    continue
                else:
                    processed[my_rsid] = 1
            	    print_bigld_output(my_rsid,rows)
         
out_matrix.close()
out_beta.close()
