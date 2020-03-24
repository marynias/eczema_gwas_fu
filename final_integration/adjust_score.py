#! /usr/bin/env python
from sys import argv
import os, sys
from collections import defaultdict as dd
import math

script, my_processed = argv

adjust_flag = 0
my_lines_dict = dd(list)
counter = 0

num_lines = len(list(open(my_processed)))
if num_lines == 1:
	sys.exit()

adj_scores = set()
with open (my_processed, 'r') as my_processed_h:
	header = my_processed_h.readline()
	my_header = header.strip().split("\t")
	first_line = my_processed_h.readline()
	counter+=1
	#Check what we are going to use: FDR, p-value, PP, score 
	my_first_line = first_line.strip().split("\t")
	#print(my_first_line)
	my_key = "line" + str(counter)
	my_lines_dict[my_key] = my_first_line
	my_index = my_header.index("FDR")
	#print(my_first_line)
	#print (my_first_line[my_index])
	if my_first_line[my_header.index("FDR")] != "NA":
		my_p_value = float(my_first_line[my_header.index("FDR")])
		if (my_p_value != 0.0 and my_p_value != 0):
			my_adjusted = math.sqrt(abs(math.log10(my_p_value)))
			adj_scores.add(my_adjusted)
		my_index = my_header.index("FDR")
		adjust_flag = 1
	elif my_first_line[my_header.index("p-value")] != "NA":
		my_p_value = float(my_first_line[my_header.index("p-value")])
		if (my_p_value != 0.0 and my_p_value != 0):
			my_adjusted = math.sqrt(abs(math.log10(my_p_value)))
			adj_scores.add(my_adjusted)
		my_index = my_header.index("p-value")
		adjust_flag = 1
	elif my_first_line[my_header.index("posterior_prob")] != "NA":
		my_p_value = my_first_line[my_header.index("posterior_prob")] 
		adj_scores.add(float(my_p_value))
		my_index = my_header.index("posterior_prob")
	elif my_first_line[my_header.index("score")] != "NA":
		my_p_value = abs(float(my_first_line[my_header.index("score")]))
		adj_scores.add(my_p_value)
		my_index = my_header.index("score")
	else:
		my_header.append("adjusted_score")
		print ("\t".join(my_header))
		my_first_line.append("1")
		print ("\t".join(my_first_line))
		for line in my_processed_h:
			lines = line.strip().split("\t")
			lines.append('1')
			print ("\t".join(lines))
		sys.exit()
	for line in my_processed_h:
		counter+=1
		lines = line.strip().split("\t")
		#print(lines)
		my_key = "line" + str(counter)
		my_lines_dict[my_key] = lines
		if lines[my_index] == 'NA':
			continue 
		my_p_value = float(lines[my_index])
		if adjust_flag:
			if (my_p_value != 0.0 and my_p_value != 0):
				my_p_value = math.sqrt(abs(math.log10(my_p_value)))
		else:
			my_p_value = abs(float(my_p_value))
		adj_scores.add(my_p_value)

#Find the biggest score.
top_value = sorted(list(adj_scores), key=float)[-1]

my_header.append("adjusted_score")
print ("\t".join(my_header))
for my_key in my_lines_dict:
	lines = my_lines_dict[my_key]
	if lines[my_index] == 'NA':
		my_increment = 1
		lines.append(str(my_increment))
		print ("\t".join(lines))
		continue
	my_p_value = float(lines[my_index])
	if adjust_flag:
		if (my_p_value != 0.0 and my_p_value != 0):
			my_p_value = math.sqrt(abs(math.log10(my_p_value)))
			my_increment = 1 + (my_p_value / float(top_value))
		else:
			my_increment = 2
	else:
		my_p_value = abs(float(my_p_value))
		my_increment = 1 + (my_p_value / float(top_value))

	lines.append(str(my_increment))
	print ("\t".join(lines))