library(biomaRt)
library(tools)

args = commandArgs(trailingOnly=TRUE)

my_results <- args[1]
my_matrix <- read.delim(my_results, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_map <- args[2]
my_map_df <- read.delim(my_map, sep="\t", stringsAsFactors=FALSE, header=TRUE)

my_matrix$se = my_matrix$beta / my_matrix$statistic
my_matrix$zscore = my_matrix$beta / my_matrix$se 

#Merge with results from Ensembl.
output <- merge(my_matrix, my_map_df,  by="gene")

out_file_gen = paste(file_path_sans_ext(my_results), "_final", sep="")
write.table(output, file=out_file_gen, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)

