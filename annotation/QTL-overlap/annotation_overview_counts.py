#! /usr/bin/env python
import os, sys, re
from sys import argv
from collections import defaultdict as dd

script, snps, my_file, output = argv
filename = os.path.basename(my_file)

def read_in_snp_list(snps):
	snp_list = {}
	with open (snps) as snps_fh:
		for line in snps_fh:
			lines = line.strip().split()
			my_snp = lines[0]	
			snp_list[my_snp] = 0
	return snp_list


def check_file(snp_list, my_file):
	MyRaw = r".*?(rs[0-9]*).*?"
	rec_pattern = re.compile(MyRaw)
	with open (my_file) as my_file_fh:
		for line in my_file_fh:
			res = rec_pattern.findall(line)
			for snp in res:
				if snp in snp_list:
					snp_list[snp] += 1
	return snp_list

out = open(output, 'w')

snp_list = read_in_snp_list(snps)
snp_list = check_file(snp_list, my_file)
all_snps = sorted(snp_list.keys())

out.write("snp" + "\t" + my_file + "\n")
for a in all_snps:
	out.write(a + "\t" + str(snp_list[a]) + "\n")
out.close()


