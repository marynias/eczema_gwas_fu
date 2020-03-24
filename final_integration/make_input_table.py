#! /usr/bin/env python
import os, sys, json, re
import argparse, gzip
from collections import defaultdict as dd

ap = argparse.ArgumentParser()

ap.add_argument('--tab',required=True,type=str,help='Input table with results of lookup.')
ap.add_argument('--tab_id',required=True,type=str,help='File table with reference id matching for input table path')
ap.add_argument('--input_path',required=True,type=str,help='Partial input path to the original input table')
ap.add_argument('--output',required=True,type=str,help='Output filename suffix')
ap.add_argument('--head',required=False,type=int,default=0,help='Line at which header starts; all previous lines ignored')
ap.add_argument('--rsid',required=True,type=int,help='Column number with rsid id in the input table')
ap.add_argument('--gene',required=True,type=int,help='Column number with gene id in the input table')
ap.add_argument('--index_locus',required=True,type=int,help='Column number with index SNP id in the input table')
ap.add_argument('--beta',required=False,type=int,default=0, help='Column number with beta values in the input table')
ap.add_argument('--beta_se',required=False,type=int,default=0, help='Column number with beta SE values in the input table')
ap.add_argument('--pvalue',required=False,type=int,default=0, help='Column number with p-value in the input table')
ap.add_argument('--fdr',required=False,type=int,default=0, help='Column number with FDR-corrected p-value in the input table')
ap.add_argument('--posterior_prob',required=False,type=int,default=0, help='Column number with posterior probability in the input table')
ap.add_argument('--bayes_factor',required=False,type=int,default=0, help='Column number with Bayes factor in the input table')
ap.add_argument('--score',required=False,type=int,default=0, help='Column number with total score in the input table')
ap.add_argument('--sig_threshold',required=False,type=float,default=0, help='Significance threshold used in the input table')
ap.add_argument('--effect_allele',required=False,type=int, help='Column number with effect allele in the input table')
ap.add_argument('--tissue',required=False,type=str,default=0, help='Header name of column with tissue/treatment in the input table, OR string with the tissue/treatment used in the input table')
ap.add_argument('--sample_size',required=False,type=str,default='NA', help='Header name of column with number of samples in the input table, OR string with the number of samples used in the input table')
ap.add_argument('--type',required=True,type=str,default='NA', help='String with the study type, e.g. eQTL, fine-mapping etc.')
ap.add_argument('--cis_trans',required=False,type=str,default=0, help='Header name of column with cis/trans status in the input table, OR string with the cis/trans used in the input table')
ap.add_argument('--evidence_weight',required=False,type=int,default=0, help='Evidence weight of data in the table. 1-prediction, 2-overlaps based on experiments, 3-statistical test based on experiments')
ap.add_argument('--delim',required=True,type=str,help='Delimiter string used in the table.')
ap.add_argument('--delim_c',required=True,type=str,help='Delimiter string used within cells with gene id.')
ap.add_argument('--delim_e',required=True,type=str,help='Delimiter string used within cells with effect_allele.')
ap.add_argument('--index_e',required=False,type=int,default=0,help='Index of the element with effect allele in the field containing it')

args = ap.parse_args()
#Initialize
my_table = args.tab
my_table_id = args.tab_id
my_output = args.output
my_rsid = args.rsid
my_gene = args.gene
head = args.head
my_index_locus = args.index_locus
my_beta = args.beta
my_beta_se = args.beta_se
my_pvalue = args.pvalue
my_pp = args.posterior_prob
my_bf = args.bayes_factor
my_score = args.score
my_sig_threshold = args.sig_threshold
my_ea = args.effect_allele
my_tissue = args.tissue
my_cis_trans = args.cis_trans
my_evidence_weight = args.evidence_weight
my_delim = args.delim
my_delim_c = args.delim_c
my_delim_e = args.delim_e
my_input_path = args.input_path
my_index_e = args.index_e
my_fdr = args.fdr
my_type = args.type
my_sample_size = args.sample_size

table_ids_dict = {}
with open (my_table_id, 'r') as my_table_id_h:
	for line in my_table_id_h:
		lines = line.strip().split()
		table_ids_dict[lines[0]] = lines[1]

#Get the table id name
out_table_id = table_ids_dict.get(my_input_path, my_input_path)
study_id = out_table_id.split("_")[0]

out = out_table_id + my_output + ".processed"
out_fh = open(out, 'w')
header = ["table_id", "study_id", "index_SNP_rsid", "current_SNP_rsid", "gene_name", "beta", "beta_se", "p-value", "FDR",
"posterior_prob", "Bayes_factor", "score", "sig_threshold", "effect_allele", "tissue", "sample_size", "study_type", "cis/trans", "evidence_weight", "sig_values_no_snp", "sig_values_no_gene"]
out_fh.write("\t".join(header) + "\n")

study_type = my_type

header_in = []
lines_to_print = dd(lambda: dd(dict))
best_pvalue = dd(lambda: dd(dict))
with open (my_table, 'r') as my_table_h:
	if head:
		for line in range(head-1):
			my_table_h.readline()
		header_in = re.split(my_delim,my_table_h.readline().strip())
	for line in my_table_h:
		lines = re.split(my_delim,line.strip())
		#print(lines)
		#print(len(lines))
		out_index_snp = lines[my_index_locus-1]
		out_rsid = lines[my_rsid-1]
		#Print seperate entry if more than one gene associated - will need to count them separately later when calculating stats for genes.
		genes = re.split(my_delim_c,lines[my_gene-1])
		if my_beta:
			out_beta = lines[my_beta-1]
		else:
			out_beta = "NA"
		if my_beta_se:
			out_beta_se = lines[my_beta_se-1]
		else:
			out_beta_se = "NA"
		if my_pvalue:
			out_pvalue = lines[my_pvalue-1]
		else:
			out_pvalue = "NA"
		if my_fdr:
			out_fdr = lines[my_fdr-1]
		else:
			out_fdr = "NA"
		if my_pp:
			out_pp = lines[my_pp-1]
		else:
			out_pp = "NA"
		if my_bf:
			out_bf = lines[my_bf-1]
		else:
			out_bf = "NA"
		if my_score:
			out_score = lines[my_score-1]
		else:
			out_score = "NA"
		if my_sig_threshold:
			out_sig_threshold = my_sig_threshold
		else:
			out_sig_threshold = "NA"
		if my_ea:
			effect = re.split(my_delim_e,lines[my_ea-1])
			out_ea = effect[my_index_e]
		else:
			out_ea = "NA"
		if my_tissue in header_in:
			my_col_t = header_in.index(my_tissue)
			out_tissue = lines[my_col_t]
		elif my_tissue:
			out_tissue = my_tissue
		else:
			out_tissue = "NA"
		if my_sample_size in header_in:
			my_col_t = header_in.index(my_sample_size)
			out_sample_size = lines[my_col_t]
		elif my_sample_size:
			out_sample_size = my_sample_size
		else:
			out_sample_size = "NA"
		if my_cis_trans:
			out_cis_trans = my_cis_trans
		else:
			out_cis_trans = "NA"
		if my_evidence_weight:
			out_evidence_weight = my_evidence_weight
		else:
			out_evidence_weight = "NA"
		for g in genes:
			my_current_line = [out_table_id, study_id, out_index_snp, out_rsid, g, str(out_beta), str(out_beta_se),
			str(out_pvalue), str(out_fdr), str(out_pp), str(out_bf), str(out_score), str(out_sig_threshold), out_ea, out_tissue, str(out_sample_size), study_type, out_cis_trans, str(out_evidence_weight)]
			out_current_line = "\t".join(my_current_line) 
			#Save only the best rsid-gene combination (lowest p-value), if possible.
			if g in lines_to_print[out_index_snp][out_rsid]:
				if out_pvalue < best_pvalue[out_index_snp][out_rsid][g]:
					best_pvalue[out_index_snp][out_rsid][g] = out_pvalue
					lines_to_print[out_index_snp][out_rsid][g] = out_current_line
			else:
				lines_to_print[out_index_snp][out_rsid][g] = out_current_line
				best_pvalue[out_index_snp][out_rsid][g] = out_pvalue
#Check how many unique significant SNP hits per locus
for index_loci in lines_to_print:
	my_unique_snps = lines_to_print[index_loci].keys()
	my_unique_genes = set()
	for my_hit in lines_to_print[index_loci]:
		my_unique_genes = set(lines_to_print[index_loci][my_hit].keys()) | my_unique_genes
	for my_hit in lines_to_print[index_loci]:
		for g in lines_to_print[index_loci][my_hit]:
			out_fh.write(lines_to_print[index_loci][my_hit][g] + "\t" + str(len(my_unique_snps)) + "\t" + str(len(my_unique_genes)) + "\n")

out_fh.close()