library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
my_significant_file <- args[2]
output <- args[3]

#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)
#Read in DEPICT result file.
my_depict <- read.delim(my_significant_file, stringsAsFactors = F, header=T)
#Filter to obtain only significantly prioritised genes.
my_depict_filtered <- my_depict[my_depict$False.discovery.rate...5. == "Yes",]
#Turn into a dataframe
length_df <- length(my_depict_filtered$Gene.symbol)
#Create a "Yes" vector
value_vector <- rep("yes", length_df)
#Combine into a df with gene names and yes as values
depict_prioritized <- data.frame(DEPICT_prioritization=value_vector, HGNC_symbol=my_depict_filtered$Gene.symbol)
#Join with the master table
combined <- merge(my_master, depict_prioritized, by="HGNC_symbol", all.x=T)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "DEPICT_prioritization")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")
