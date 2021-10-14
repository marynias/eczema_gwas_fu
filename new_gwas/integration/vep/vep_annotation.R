library("tools")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
my_significant_file <- args[2]
output <- args[3]

my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)
#Read in DEPICT result file.
my_vep <- read.delim(my_significant_file, stringsAsFactors = F, header=T)
#Filter to contain only intron-containing variants
my_vep_intron <- my_vep[my_vep$Consequence == "intron_variant",]
#Filter to contain only missense variants
my_vep_missense <- my_vep[my_vep$Consequence == "missense_variant",]
#Turn into a dataframe
length_intron <- length(unique(my_vep_intron$Gene))
#Create a "Yes" vector
value_vector_intron <- rep("yes", length_intron)

#Turn into a dataframe
length_missense <- length(unique(my_vep_missense$Gene))
#Create a "Yes" vector
value_vector_missense <- rep("yes", length_missense)

#Combine into a df with gene names and yes as values
intron_prioritized <- data.frame(VEP_intron=value_vector_intron, Ensembl_gene_ID=unique(my_vep_intron$Gene))
missense_prioritized <- data.frame(VEP_missense=value_vector_missense, Ensembl_gene_ID=unique(my_vep_missense$Gene))

#Join with the master table
combined <- merge(my_master, intron_prioritized, by="Ensembl_gene_ID", all.x=T)
combined2 <- merge(combined, missense_prioritized, by="Ensembl_gene_ID", all.x=T)

#Resort the table.
combined2 <- combined2[gtools::mixedorder(combined2$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "VEP_intron", "VEP_missense")
combined2 <- combined2[my_col_order]
#Write table to file
write.table(combined2, output, quote=F, sep=",", row.names=F, na="")