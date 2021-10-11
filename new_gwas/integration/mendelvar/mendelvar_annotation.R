library("tools")
library("stringr")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
my_significant_file <- args[2]
output <- args[3]

#Read in all results files
temp = list.files(pattern="*parsed")
#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)
my_files = lapply(temp, read.delim, stringsAsFactors = F,)
final_df = do.call(rbind, my_files)
#Find significant pathways
significant_df <- final_df[final_df$Bootstrapped_p_value < 0.06,]
#Save the results to file
write.table(significant_df, my_significant_file, quote=F, sep="\t", row.names=F)
#Obtain the names of the genes with significant enrichment
list_genes <- significant_df$Overlapping_genes_hgnc_names
listed_genes <- as.list(strsplit(list_genes, ","))
listed_genes <- unlist(listed_genes)
#Remove whitespace
listed_genes <- gsub("\\s", "", listed_genes) 
#Remove duplicate names
final_list <- unique(listed_genes)
#Turn into a dataframe
length_df <- length(final_list)
#Create a "Yes" vector
value_vector <- rep("yes", length_df)
#Combine into a df with gene names and yes as values
mendelvar_enrichment <- data.frame(MendelVar_sig_enrichment=value_vector, HGNC_symbol=final_list)
#Join with the master table
combined <- merge(my_master, mendelvar_enrichment, by="HGNC_symbol", all.x=T)
#Find terms with "kera", "derma", "skin" keywords among all the DO and HPO results.
#Get HPO and DO files. 
my_input <- c("hpo.out.inrich.parsed", "do.out.inrich.parsed")
my_in = lapply(my_input, read.delim, stringsAsFactors = F,)
do_hpo = do.call(rbind, my_in)
my_search_terms <- c("skin", "kera", "derma")
skin_matches <- str_detect(do_hpo$Gene_set_description, regex(paste(my_search_terms, collapse = '|'), ignore_case = TRUE))
skin_matches_df <- do_hpo[skin_matches,]
#Obtain the names of the genes with the keywords
list_genes <- skin_matches_df$Overlapping_genes_hgnc_names
listed_genes <- as.list(strsplit(list_genes, ","))
listed_genes <- unlist(listed_genes)
#Remove whitespace
listed_genes <- gsub("\\s", "", listed_genes) 
#Remove duplicate names
final_list <- unique(listed_genes)
#Turn into a dataframe
length_df <- length(final_list)
#Create a "Yes" vector
value_vector <- rep("yes", length_df)
#Combine into a df with gene names and yes as values
mendelvar_keywords <- data.frame(MendelVar_skin_keywords=value_vector, HGNC_symbol=final_list)
#Join with the master table
combined2 <- merge(combined, mendelvar_keywords, by="HGNC_symbol", all.x=T)
#Resort the table.
combined2 <- combined2[gtools::mixedorder(combined2$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "MendelVar_sig_enrichment",
                  "MendelVar_skin_keywords")
combined <- combined[my_col_order]
#Write table to file
write.table(combined2, "paternoster2015_mendelvar.csv", quote=F, sep=",", row.names=F, na="")



