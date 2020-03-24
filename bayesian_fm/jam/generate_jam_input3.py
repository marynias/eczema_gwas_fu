#! /usr/bin/env python
import os, sys, re, argparse
import gzip
from collections import defaultdict as dd
ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with GWAS summary stats')
ap.add_argument('--vcf',required=True,type=str,help='Input file with reference file genotypes')
ap.add_argument('--out',required=True,type=str,help='Prefix for JAM input files')
ap.add_argument('--pos',required=True,type=int,help='Column with SNP pos in --tab')
ap.add_argument('--rsid',required=True,type=int,help='Column with SNP rsid in --tab')
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
vcf = args.vcf
pos = args.pos
ref = args.ref
alt = args.alt
beta = args.beta
se = args.se
rsid = args.rsid


#Store effect and non-effect allele identity, as well as the calculated z-score.
rsids = dd(lambda: dd(str))
with open (tab) as tab_h:
    headers = tab_h.readline()
    for line in tab_h:
        lines = line.strip().split()  
        rsids[lines[rsid-1]]['ref'] = lines[ref-1]
        rsids[lines[rsid-1]]['alt'] = lines[alt-1]
        rsids[lines[rsid-1]]['beta'] = lines[beta-1]
        rsids[lines[rsid-1]]['se'] = lines[se-1]

out_matrix = open(out + ".matrix", 'w')
out_beta = open(out + ".beta", 'w')
out_beta.write("rsid" + "\t" + "beta" + "\t" + "se" + "\n")
processed = rsids.keys()

def print_jam_output(my_rsid, rows, switch):
    out_beta.write(my_rsid + "\t" + str(rsids[my_rsid]['beta']) + "\t" + str(rsids[my_rsid]['se']) + "\n")
    out_matrix.write(my_rsid)
    if switch == 1:
        count_val = "1"
    else:
        count_val = "0"
    for indin in rows[9:]:
        haplos = indin.split(":")
        if haplos[0] != "./.":
            genotypes = haplos[0].split("/")
            my_count = genotypes.count(count_val)
            out_matrix.write("\t" + str(my_count))
        else:
             out_matrix.write("\t" + "NA")
    out_matrix.write("\n")

#Iterate through VCF files to print genotype matrix and beta vector.
file_stream = gzip.open(vcf, 'rt')
for lines in file_stream:
    if(lines[0:2] != "##"):
        if (lines[0:6] == "#CHROM"):
            header = lines.strip().split()
            A0_index = header.index("REF")
            A1_index = header.index("ALT")
            rsid_index = header.index("ID")
            out_matrix.write("rsid" + "\t" + "\t".join(header[9:]) + "\n")
        else:
            rows = lines.strip().split()
            my_rsids = rows[rsid_index].split(",")
            my_rsid = my_rsids[0]
            if str(my_rsid) in processed:
                vcf_A0 = rows[A0_index]
                vcf_A1 = rows[A1_index]
                if vcf_A1 == rsids[my_rsid]['alt'] and vcf_A0 == rsids[my_rsid]['ref']:
                    print_jam_output(my_rsid,rows,0)
                elif vcf_A1 == rsids[my_rsid]['ref'] and vcf_A0 == rsids[my_rsid]['alt']:
                    print_jam_output(my_rsid,rows,1)
                else:
                    continue
out_beta.close()
out_matrix.close()
file_stream.close()