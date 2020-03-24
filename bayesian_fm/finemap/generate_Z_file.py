#! /usr/bin/env python
import os, sys, re, argparse
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with GWAS summary stats')
ap.add_argument('--proces',required=True,type=str,help="Table with SNPs processed to generate LD matrix with Paintor tool")
ap.add_argument('--out',required=True,type=str,help='Name of the output loc file')
ap.add_argument('--chrom',required=True,type=int,help='Column with chromosome number')
ap.add_argument('--pos',required=True,type=int,help='Column with SNP position')
ap.add_argument('--ident',required=True,type=int,help='Column with SNP id')
ap.add_argument('--ref',required=True,type=int,help='Column with ref allele')
ap.add_argument('--alt',required=True,type=int,help='Column with alf allele')
ap.add_argument('--beta',required=True,type=int,help='Column with effect size')
ap.add_argument('--se',required=True,type=int,help='Column with SE')
ap.add_argument('--eas',required=True,type=int,help='Column with effect allele frequency')

args = ap.parse_args()
#Initialize
tab = args.tab
out = args.out
proces = args.proces
chrom = args.chrom
pos = args.pos
ident = args.ident
ref = args.ref
alt = args.alt
beta = args.beta
se = args.se 
eas = args.eas

#Extract SNP ids to be used in downstream analysis:
to_keep = {}
proces_h = open(proces)
for line in proces_h:
    lines = line.split()
    to_keep[lines[2]] = 1
proces_h.close()

tab_h = open(tab)
out_h = open(out, 'w')
out_h.write("rsid chromosome position allele1 allele2 maf beta se" + "\n")
headers = tab_h.readline()
for line in tab_h:
    lines = line.split()    
    my_ident = lines[ident-1]
    if my_ident in to_keep:
        my_chrom = int(lines[chrom-1])
        my_pos = lines[pos-1]
        my_ref = lines[ref-1]
        my_alt = lines[alt-1]
        my_beta = lines[beta-1]
        my_se = lines[se-1]
        my_eas = float(lines[eas-1])
        if my_eas > 0.5:
            my_maf = 1 - my_eas
        else:
            my_maf = my_eas
        out_h.write(my_ident + " " + str(my_chrom) + " " + str(my_pos) + " ",)
        out_h.write(my_ref + " " + my_alt + " " + str(my_maf) + " ",)
        out_h.write(str(my_beta) + " " + str(my_se) + "\n")
    else:
        continue
tab_h.close()
out_h.close()