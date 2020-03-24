#! /usr/bin/env python
import os, sys, re, argparse
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with GWAS summary stats')
ap.add_argument('--out',required=True,type=str,help='Name of the output loc file')
ap.add_argument('--chrom',required=True,type=int,help='Column with chromosome number')
ap.add_argument('--chr_out',required=True,type=int,help='Chromosome no to output')
ap.add_argument('--pos',required=True,type=int,help='Column with SNP position')
ap.add_argument('--ident',required=True,type=int,help='Column with SNP id')
ap.add_argument('--ref',required=True,type=int,help='Column with ref allele')
ap.add_argument('--alt',required=True,type=int,help='Column with alf allele')
ap.add_argument('--beta',required=True,type=int,help='Column with effect size')
ap.add_argument('--se',required=True,type=int,help='Column with SE')
ap.add_argument('--indels',required=True,type=str,help='Output indels? Y/N', default='N')

args = ap.parse_args()

#Initialize
tab = args.tab
out = args.out
chrom = args.chrom
chr_out = args.chr_out
pos = args.pos
ident = args.ident
ref = args.ref
alt = args.alt
beta = args.beta
se = args.se 
indels = args.indels

tab_h = open(tab)
out_h = open(out, 'w')
out_h.write("chr pos rsid A0 A1 Zscore"  + "\n")
headers = tab_h.readline()
for line in tab_h:
    lines = line.split()
    #print(lines)
    my_chrom = int(lines[chrom-1])
    my_pos = lines[pos-1]
    my_ident = lines[ident-1]
    my_ref = lines[ref-1]
    my_alt = lines[alt-1]
    my_beta = lines[beta-1]
    my_se = lines[se-1]
    if (not my_ident.startswith("rs") and indels == "N"):
        continue
    if my_chrom == chr_out:
        Zscore = float(my_beta) / float(my_se)
        out_h.write(str(my_chrom) + " " + str(my_pos) + " " + my_ident + " ",)
        out_h.write(my_ref + " " + my_alt + " " + str(Zscore) + "\n")        
tab_h.close()
out_h.close()