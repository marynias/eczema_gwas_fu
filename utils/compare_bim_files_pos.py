#! /usr/bin/env python

#Compare two files by position
#Allele values have to be the same case (upper or lower) in both files being compared.

import os, sys
from sys import argv
from collections import defaultdict as dd
from Bio.Seq import Seq

script, file1, file2 = argv
filename1 = os.path.basename(file1)
filename2 = os.path.basename(file2)

loci = dd(lambda: dd(list))
rsid = {}
loci2 = dd(lambda: dd(list))
rsid2 = {}

with open (file1) as fh1:
	for line in fh1:
		lines = line.strip().split()
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		loci[chrom][pos] = [snp_id, allele_1, allele_2]
		rsid[snp_id] = [chrom, pos, allele_1, allele_2]

with open (file2) as fh2:
	for line in fh2:
		lines = line.strip().split()
		#print(lines)
		chrom = lines[0]
		pos = lines[3]
		snp_id = lines[1]
		allele_1 = lines[4]
		allele_2 = lines[5]
		loci2[chrom][pos] = [snp_id, allele_1, allele_2]
		rsid2[snp_id] = [chrom, pos, allele_1, allele_2]

#Check how many rsids present in file 1, present in file 2
def rsid_check(rsid, rsid2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_missing"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present"
	out1_fh = open (out1, 'w')
	out2_fh = open (out2, 'w')
	missing = set()
	present = set()
	positives = 0
	counter = 0
	for r in rsid:
		counter += 1
		if r in rsid2:
			positives += 1
			present.add(r)
			out2_fh.write(r + "\n")
		else:
			missing.add(r)
			out1_fh.write(r + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	print ("Total number of rsid IDs found in %s matching in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants in the file1, that is %.2f percent." % my_percentage + " --> rsid_check")
	out1_fh.close()
	out2_fh.close()
	return (missing, present)

#Out of rsids found in both files, check how many share the same position
def pos_check(present_rsids, rsid, rsid2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_present.pos_diff"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present.pos_match"
	out1_fh = open(out1, 'w')
	out2_fh = open(out2, 'w')
	differences = set()
	matching = set()
	positives = 0
	counter = 0
	for r in present_rsids:
		counter += 1
		pos1 = rsid[r][1]
		pos2 = rsid2[r][1]
		if pos1 == pos2:
			positives += 1
			matching.add(r)
			out2_fh.write(r + "\n")
		else:
			differences.add(r)
			out1_fh.write(r + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	print ("Total number of variant positions found in %s matching in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with the same rsid, that is %.2f percent." % my_percentage + " --> pos_check")
	out1_fh.close()
	out2_fh.close()
	return (differences, matching)

#Out of rsids found in both files and matching positions, check how many share the same alleles (order not important), and how many of them are ambiguous?
def allele_check_rsid (matching_pos_rsid, rsid, rsid2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_present.allele_diff"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present.allele_match"
	out3 = filename1 + "_vs_" + filename2 + ".rsid_present.allele_match.ambig"
	out1_fh = open(out1, 'w')
	out2_fh = open(out2, 'w')
	out3_fh = open(out3, 'w')
	differences = set()
	matching = set()
	positives = 0
	ambig = 0
	counter = 0
	for r in matching_pos_rsid:
		counter += 1
		allele1_1 = rsid[r][2]
		allele1_2 = rsid[r][3]
		allele2_1 = rsid2[r][2]
		allele2_2 = rsid2[r][3]
		if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
			positives += 1
			matching.add(r)
			out2_fh.write(r + "\n")
			if allele1_1 == str(Seq(allele2_2).reverse_complement()) and allele1_2 == str(Seq(allele2_1).reverse_complement()):
				ambig += 1
				out3_fh.write(r + "\n")
		elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
			positives += 1
			matching.add(r)
			out2_fh.write(r + "\n")
			if allele1_1 == str(Seq(allele2_1).reverse_complement()) and allele1_2 == str(Seq(allele2_2).reverse_complement()):
				ambig += 1
				out3_fh.write(r + "\n")
		else:
			differences.add(r)
			out1_fh.write(r + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	if positives > 0:
		ambig_percentage = (float(ambig) / float(positives)) * 100
	else:
		ambig_percentage = 0
	print ("Total number of variant alleles found in %s matching in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with the same rsid and position, that is %.2f percent. " % my_percentage + "{:,}".format(ambig) + " out of " "{:,}".format(positives) + " variants matched for alleles are ambiguous, that is %.2f percent." % ambig_percentage + " --> allele_check_rsid")
	out1_fh.close()
	out2_fh.close()
	out3_fh.close()
	return (differences, matching)

#Check if position occupied by rsids which were not found in file2, as per routine rsid_check, are occupied by 
#variants with different snp_id in file2.
def missing_rsid_check(missing_rsids, rsid, loci2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_missing.pos_nomatch"
	out1_fh = open(out1, 'w')
	out2 = filename1 + "_vs_" + filename2 + ".rsid_missing.pos_match"
	out2_fh = open(out2, 'w')
	out3 = filename1 + "_vs_" + filename2 + ".rsid_missing.pos_nomatch_chrom"
	out3_fh = open(out3, 'w')
	positives = 0
	counter = 0
	matching = set()
	differences = set()
	for r in missing_rsids:
		counter += 1
		my_chrom = rsid[r][0]
		my_pos = rsid[r][1]
		if my_chrom in loci2:
			if my_pos in loci2[my_chrom]:
				positives += 1
				other_rsid = (my_chrom, my_pos)
				matching.add(other_rsid)
				out2_fh.write(my_chrom + "\t" + my_pos + "\n")
			else:
				other_rsid = (my_chrom, my_pos)
				differences.add(other_rsid)
				out1_fh.write(r + "\n")
		else:
			other_rsid = (my_chrom, my_pos)
			differences.add(other_rsid)
			out3_fh.write(r + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	print ("Total number of alternative rsids found in %s matching for position in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with different rsid, that is %.2f percent." % my_percentage + " --> missing_rsid_check")
	out1_fh.close()
	out2_fh.close()
	out3_fh.close()
	return (differences, matching)

#Check for matching alleles when comparing by variant position. See output of subroutine above.
def allele_check_loci (matching_pos_different_rsids, loci, loci2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_missing.pos_match.allele_diff"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present.pos.match.allele_match"
	out3 = filename1 + "_vs_" + filename2 + ".rsid_present.pos.match.allele_match.ambig"
	out1_fh = open(out1, 'w')
	out2_fh = open(out2, 'w')
	out3_fh = open(out3, 'w')
	positives = 0
	counter = 0
	ambig = 0
	positions = set()
	differences = set()
	for r in matching_pos_different_rsids:
		my_chrom = r[0]
		my_pos = r[1]
		allele1_1 = loci[my_chrom][my_pos][1]
		allele1_2 = loci[my_chrom][my_pos][2]
		allele2_1 = loci2[my_chrom][my_pos][1]
		allele2_2 = loci2[my_chrom][my_pos][2]
		counter += 1
		if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
			positives += 1
			rsids = (loci[my_chrom][my_pos][0],loci2[my_chrom][my_pos][0])
			positions.add(rsids)
			out2_fh.write(loci[my_chrom][my_pos][0] + "\t" + loci2[my_chrom][my_pos][0] + "\n")
			if allele1_1 == str(Seq(allele2_2).reverse_complement()) and allele1_2 == str(Seq(allele2_1).reverse_complement()):
				ambig += 1
				out3_fh.write(loci[my_chrom][my_pos][0] + "\t" + loci2[my_chrom][my_pos][0] + "\n")
		elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
			positives += 1
			rsids = (loci[my_chrom][my_pos][0],loci2[my_chrom][my_pos][0])
			positions.add(rsids)
			out2_fh.write(loci[my_chrom][my_pos][0] + "\t" + loci2[my_chrom][my_pos][0] + "\n")
			if allele1_1 == str(Seq(allele2_1).reverse_complement()) and allele1_2 == str(Seq(allele2_2).reverse_complement()):
				ambig += 1
				out3_fh.write(loci[my_chrom][my_pos][0] + "\t" + loci2[my_chrom][my_pos][0] + "\n")
		else:
			rsids = (loci[my_chrom][my_pos][0],loci2[my_chrom][my_pos][0])
			differences.add(rsids)
			out1_fh.write(loci[my_chrom][my_pos][0] + "\t" + loci2[my_chrom][my_pos][0] + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	if positives > 0:
		ambig_percentage = (float(ambig) / float(positives)) * 100
	else:
		ambig_percentage = 0
	print ("Total number of alternative rsids found in %s matching for alleles in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with the same position and different rsid, that is %.2f percent. " % my_percentage + "{:,}".format(ambig) + " out of " "{:,}".format(positives) + " variants matched for alleles are ambiguous, that is %.2f percent." % ambig_percentage + " --> allele_check_loci")
	out1_fh.close()
	out2_fh.close()
	out3_fh.close()
	return (differences, positions)
#Check if any alleles not matching are potentially due to flips in strand used, for snps with matching rsids.
def flip_check_rsid (different_allele_same_rsid_same_pos, rsid, rsid2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_present.allele_diff.noallele_flip"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present.allele_diff.allele_flip"
	out1_fh = open(out1, 'w')
	out2_fh = open(out2, 'w')
	differences = set()
	matching = set()
	positives = 0
	counter = 0
	for r in different_allele_same_rsid_same_pos:
		counter += 1
		allele1_1 = rsid[r][2]
		allele1_2 = rsid[r][3]
		#Flip alleles in the second file.
		allele2_1 = str(Seq(rsid2[r][2]).reverse_complement())
		allele2_2 = str(Seq(rsid2[r][3]).reverse_complement())
		if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
			positives += 1
			matching.add(r)
			out2_fh.write(r + "\n")
		elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
			positives += 1
			matching.add(r)
			out2_fh.write(r + "\n")
		else:
			differences.add(r)
			out1_fh.write(r + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	print ("Total number of rsid IDs found in %s potentially matching when alleles flipped in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with the same rsid and position, that is %.2f percent. " % my_percentage + " --> flip_check_rsid")
	out1_fh.close()
	out2_fh.close()
	return (differences, matching)

#Check if any alleles not matching are potentially due to flips in strand used, for snps with not matching rsids, ut matching chrom and position.
def flip_check_pos(different_allele_different_rsid_same_pos, rsid, rsid2):
	out1 = filename1 + "_vs_" + filename2 + ".rsid_missing.pos_match.allele_diff.noallele_flip"
	out2 = filename1 + "_vs_" + filename2 + ".rsid_present.pos.match.allele_match.allele_flip"
	out1_fh = open(out1, 'w')
	out2_fh = open(out2, 'w')
	positives = 0
	counter = 0
	matching = set()
	differences = set()
	for my_set in different_allele_different_rsid_same_pos:
		r1 = my_set[0]
		r2 = my_set[1]
		counter += 1
		allele1_1 = rsid[r1][2]
		allele1_2 = rsid[r1][3]
		#Flip alleles in the second file.
		allele2_1 = str(Seq(rsid2[r2][2]).reverse_complement())
		allele2_2 = str(Seq(rsid2[r2][3]).reverse_complement())
		if (allele1_1 == allele2_1 and allele1_2 == allele2_2):
			positives += 1
			rsids = (r1,r2) 
			matching.add(rsids)
			out2_fh.write(r1 + "\t" + r2 + "\n")

		elif (allele1_2 == allele2_1 and allele1_1 == allele2_2):
			positives += 1
			rsids = (r1,r2) 
			matching.add(rsids)
			out2_fh.write(r1 + "\t" + r2 + "\n")
		else:
			rsids = (r1,r2) 
			differences.add(rsids)
			out1_fh.write(r1 + "\t" + r2 + "\n")
	if counter > 0:
		my_percentage = (float(positives) / counter) * 100
	else:
		my_percentage = 0
	print ("Total number of alternative rsids found in %s potentially matching for alleles when alleles flipped in %s is " % (filename1, filename2) + "{:,}".format(positives) + " out of " "{:,}".format(counter) + " variants with the same position and different rsid, that is %.2f percent. " % my_percentage + " --> flip_check_pos")
	out1_fh.close()
	out2_fh.close()
	return (differences, matching)

missing_rsids, present_rsids = rsid_check(rsid, rsid2)
different_pos_same_rsid, matching_pos_same_rsid = pos_check(present_rsids, rsid, rsid2)
different_allele_same_rsid_same_pos, matching_allele_same_rsid_same_pos = allele_check_rsid(matching_pos_same_rsid, rsid, rsid2)
if len(different_allele_same_rsid_same_pos) > 0:
	non_flipped_different_allele_same_rsid_same_pos, flipped_different_allele_same_rsid_same_pos = flip_check_rsid(different_allele_same_rsid_same_pos, rsid, rsid2)

if len(missing_rsids) > 0:
	nonmatching_pos_different_rsids, matching_pos_different_rsids = missing_rsid_check(missing_rsids, rsid, loci2)
	if len(matching_pos_different_rsids) > 0:
		different_allele_different_rsid_same_pos, matching_allele_different_rsid_same_pos = allele_check_loci(matching_pos_different_rsids, loci, loci2)
		if len(different_allele_different_rsid_same_pos) > 0:
			non_flipped_different_allele_different_rsid_same_pos, flipped_different_allele_different_rsid_same_pos = flip_check_pos(different_allele_different_rsid_same_pos, rsid, rsid2)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             