#! /usr/bin/env python

import os, sys, re, argparse, gzip

#Check if the Ensembl transcript to probe mapping in files 
#ProbeAnnotation2_CorrectlyMapped_mapping_annotations_added_20170301.txt (the main reference)
#and eQTLsFDR-ProbeLevel.txt.gz matches.
#Can be used to compare other types of files with ID values matched with another value.

ap = argparse.ArgumentParser()
ap.add_argument('--tab',required=True,type=str,help='Input table with eQTL summary stats')
ap.add_argument('--proces',required=True,type=str,help='Mapping file')
ap.add_argument('--tab_keys',required=True,type=int,help='Column with IDs in the tab file')   
ap.add_argument('--tab_vals',required=True,type=int,help='Column with values in the tab file')  
ap.add_argument('--proces_keys',required=True,type=int,help='Column with IDs in the proces file') 
ap.add_argument('--proces_vals',required=True,type=int,help='Column with values in the proces file')   

args = ap.parse_args()
#Initialize
tab = args.tab
proces = args.proces
tab_keys = args.tab_keys
tab_vals = args.tab_vals
proces_keys = args.proces_keys
proces_vals = args.proces_vals


loci = {}
with open (proces) as proces_h:
	header = proces_h.readline()
	for line in proces_h:
		lines = line.strip().split("\t")
		loci[lines[proces_keys-1]] = [a.strip() for a in lines[proces_vals-1].split(";")].sort()

if tab.endswith(".gz"):
	tab_h = gzip.open(tab, 'rt')
else:
	tab_h = open(tab, 'r')

header = tab_h.readline()
for line in tab_h:
	lines = line.strip().split("\t")
	if lines[tab_keys-1] in loci:
		current_vals = [a.strip() for a in lines[tab_vals-1].split(",")].sort()
		if current_vals == loci[lines[tab_keys-1]]:
			pass
		else:
			print ("Warning, different values for key %s: %s in --tab and %s in --proces" % (lines[tab_keys-1], ",".join(current_vals), ",".join(loci[lines[tab_keys-1]])))
	else:
		print("Warning, key %s not present in the mapping file" % lines[tab_keys-1])

tab_h.close()

