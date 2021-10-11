library("dplyr")
library("stringr")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
my_significant_file <- args[2]
output <- args[3]


#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

#Read in MAGMA gene output
my_magma <- read.delim(my_significant_file, stringsAsFactors = F, header=T, sep="")

#Get the genes below alpha value with Bonferroni correction
alpha <- 0.05/dim(my_magma)[1]

my_magma <- as_tibble(my_magma) %>% filter(P < alpha)
#Turn into a dataframe
length_df <- length(my_magma$GENE)
#Create a "Yes" vector
value_vector <- rep("yes", length_df)
#Combine into a df with gene names and yes as values
magma_prioritized <- data.frame(MAGMA_prioritization=value_vector, Ensembl_gene_ID=my_magma$GENE)
#Join with the master table
combined <- merge(my_master, magma_prioritized, by="Ensembl_gene_ID", all.x=T)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "MAGMA_prioritization")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")
