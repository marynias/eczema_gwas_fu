#! /usr/bin/env python
import os, sys, json, re
import argparse, gzip
from collections import defaultdict as dd

ap = argparse.ArgumentParser()

ap.add_argument('--tab',required=False,default=sys.stdin,type=argparse.FileType('r'), nargs='?', help='Input table with results. Can be streamed to script. ')
ap.add_argument('--db',required=True,type=str,help='Input JSON dictionary matching aliases to dbSNP reference database / HUGO gene ref database')
ap.add_argument('--head',required=False,type=int,default=0,help='Line at which header starts; all previous lines ignored')
ap.add_argument('--my_id',required=True,type=int,help='Column number with rsid id / gene_id in the input table')
ap.add_argument('--delim',required=True,type=str,help='Delimiter string used in the table.')
ap.add_argument('--delim_c',required=True,type=str,help='Delimiter string used within cells.')
ap.add_argument('--convert',required=False,type=str,help='Convert the rsid id/gene id to lowercase or uppercase? - arguments lowercase, uppercase accepted')

args = ap.parse_args()
#Initialize
fh2 = args.tab
db = args.db
my_rsid = args.my_id
head = args.head
my_delim = args.delim
my_delim_c = args.delim_c
my_convert = args.convert

if my_convert:
	if my_convert == "upper":
		def convert_ids(values):
			converted = [db_dict.get(v.strip().upper(), v.strip().upper()) for v in values]
			return converted
	elif my_convert == "lower":
		def convert_ids(values):
			converted = [db_dict.get(v.strip().lower(), v.strip().lower()) for v in values]
			return converted
else:
	def convert_ids(values):
		converted = [db_dict.get(v.strip(), v.strip()) for v in values]
		return converted


#Load in dbSNP/gene_names mappings.
with open (db, 'r') as db_h:
	db_dict = json.load(db_h)

if head:
	for line in range(head-1):
		fh2.readline()
	header = fh2.readline()
	print (header.strip())
for line in fh2:
	lines = re.split(my_delim,line.strip())
	#print (lines)
	if len(lines) >= my_rsid:
		values = re.split(my_delim_c,lines[my_rsid-1])
		#print (values)
		substituted = convert_ids(values)
		out_v = my_delim_c.join(substituted)
		lines[my_rsid-1] = out_v
		print (my_delim.join(lines))
	else:
		pass
