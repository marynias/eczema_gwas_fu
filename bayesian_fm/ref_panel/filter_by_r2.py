#! /usr/bin/env python
import os, sys, re
from collections import defaultdict as dd
from sys import argv
import gzip

#Filter VCF output with LD calculation to output only positions with r2 higher or equal to X relative to SNP Y.
script, my_pos, my_r2, my_ld_file = argv

if my_ld_file.endswith(".gz"):
	file_stream = gzip.open(my_ld_file, 'rt')
else:
	file_stream = open(my_ld_file, 'r')
for line in file_stream:
    lines = line.strip().split()
    if (lines[1] == my_pos and float(lines[4]) >= float(my_r2)):
    	print (lines[2], "\t", lines[4])
    elif (lines[2] == my_pos and float(lines[4]) >= float(my_r2)):
    	print (lines[1], "\t", lines[4])
