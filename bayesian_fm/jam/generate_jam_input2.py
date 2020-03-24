#! /usr/bin/env python
import os, sys, re, argparse
import gzip
from collections import defaultdict as dd
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with GWAS summary stats')
ap.add_argument('--proces',required=True,type=str,help="Table with SNPs processed to generate LD matrix for finemap")
ap.add_argument('--vcf',required=True,type=str,help='Input file with reference file genotypes')
ap.add_argument('--out',required=True,type=str,help='Prefix for JAM input files')
ap.add_argument('--pos',required=True,type=int,help='Column with SNP pos in --tab')
ap.add_argument('--ref',required=True,type=int,help='Column with ref allele in --tab')
ap.add_argument('--alt',required=True,type=int,help='Column with alf allele in --tab')
ap.add_argument('--se',required=True,type=int,help='Column with SE in --tab')
ap.add_argument('--beta',required=True,type=int,help='Column with effect size in --tab')

#The script assumes all missing data has been imputed!
#The script assumes independent file for each chrom.
#Outputting beta and se, rather than z-score for each processed SNP.

args = ap.parse_args()
#Initialize
tab = args.tab
out = args.out
proces = args.proces
vcf = args.vcf
pos = args.pos
ref = args.ref
alt = args.alt
beta = args.beta
se = args.se

#Dictionary to store positions to be processed
rsids = dd(lambda: dd(str))
with open (proces) as proces_h:
    headers = proces_h.readline()
    for line in proces_h:
        lines = line.strip().split()  
        rsids[lines[1]]['proc'] = lines[2]

#Store effect and non-effect allele identity, as well as the calculated z-score.
with open (tab) as tab_h:
    headers = tab_h.readline()
    for line in tab_h:
        lines = line.strip().split()  
        if lines[pos-1] in rsids:
            rsids[lines[pos-1]]['ref'] = lines[ref-1]
            rsids[lines[pos-1]]['alt'] = lines[alt-1]
            rsids[lines[pos-1]]['beta'] = lines[beta-1]
            rsids[lines[pos-1]]['se'] = lines[se-1]
        else:
            continue
out_matrix = open(out + ".matrix", 'w')
out_beta = open(out + ".beta", 'w')
out_beta.write("rsid" + "\t" + "beta" + "\t" + "se" + "\n")
processed = rsids.keys()

def print_jam_output(my_pos, rows, switch):
    out_beta.write(rsids[my_pos]['proc'] + "\t" + str(rsids[my_pos]['beta']) + "\t" + str(rsids[my_pos]['se']) + "\n")
    out_matrix.write(rsids[my_pos]['proc'])
    if switch == 1:
        count_val = "1"
    else:
        count_val = "0"
    for haplos in rows[9:]:
        genotypes = haplos.split("|")
        my_count = genotypes.count(count_val)
        out_matrix.write("\t" + str(my_count))
    out_matrix.write("\n")

#Iterate through VCF files to print genotype matrix and beta vector.
file_stream = gzip.open(vcf, 'rt')
for lines in file_stream:
    if(lines[0:2] != "##"):
        if (lines[0:6] == "#CHROM"):
            header = lines.strip().split()
            A0_index = header.index("REF")
            A1_index = header.index("ALT")
            pos_index = header.index("POS")
            out_matrix.write("rsid" + "\t" + "\t".join(header[9:]) + "\n")
        else:
            rows = lines.strip().split()
            my_pos = rows[pos_index]
            if str(my_pos) in processed:
                vcf_A0 = rows[A0_index]
                vcf_A1 = rows[A1_index]
                if vcf_A1 == rsids[my_pos]['alt'] and vcf_A0 == rsids[my_pos]['ref']:
                    print_jam_output(my_pos,rows,0)
                elif vcf_A1 == rsids[my_pos]['ref'] and vcf_A0 == rsids[my_pos]['alt']:
                    print_jam_output(my_pos,rows,1)
                else:
                    continue
out_beta.close()
out_matrix.close()
file_stream.close()