#! /usr/bin/env python
import os, sys, re, argparse
from collections import defaultdict as dd
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with GWAS summary stats')
ap.add_argument('--ped',required=True,type=str,help="Transposed ped file produced by Plink")
ap.add_argument('--out',required=True,type=str,help='Name of the output file')
ap.add_argument('--chrom',required=True,type=int,help='Column with chromosome number in --tab file')
ap.add_argument('--pos',required=True,type=int,help='Column with SNP position in --tab file')
ap.add_argument('--ident',required=True,type=int,help='Column with SNP id in --tab file')

args = ap.parse_args()

snptable = args.tab
ped = args.ped
out_f = args.out
chrom = args.chrom
pos = args.pos
ident = args.ident

#Dictionary to store positions to be processed
rsids = dd(lambda: dd(str))
with open (snptable) as snptable_h:
    headers = snptable_h.readline()
    for line in snptable_h:
        lines = line.strip().split()
        my_ident = lines[ident-1]
        my_chrom = lines[chrom-1]
        my_pos = lines[pos-1] 
        rsids[my_chrom][my_pos] = my_ident

out = open(out_f, 'w')
 #Iterate over the input tped file.
with open (ped) as ped_h:
	for line in ped_h:
		lines = line.strip().split()  
		if rsids[lines[0]][lines[3]] != "": 
			#Check that positions match (not checking for alleles, since this is done in the intial
			#preprocessing step by Paintor's script sub_CalcLD_1KG_VC_1K.sh 
			if rsids[lines[0]][lines[3]] != lines[1]:
				print ("Warning, %s in the tped file is under %s in the GWAS file" % (lines[1], rsids[lines[0]][lines[3]]))
				lines[1] = rsids[lines[0]][lines[3]]
				out.write(" ".join(lines) + "\n")
			else:
				out.write(line)
out.close()
