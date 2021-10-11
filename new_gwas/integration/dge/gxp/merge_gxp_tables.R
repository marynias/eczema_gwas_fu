library("stringr")
library("dplyr")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
output <- args[2]

#Merge GXP tables to master table.

my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

#Read in GXP tables and merge them
temp = list.files(pattern="*_gxp.tsv")

loadFile <- function(x) {
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  df
}

all_normal <- lapply(temp, loadFile)
final_df = as_tibble(do.call(rbind, all_normal))

#Temprorary join of the two files by ENSG
temp_join <- merge(my_master, final_df, by="Ensembl_gene_ID")
#Temprorary join of the two files by symbol
temp_join2 <- merge(my_master, final_df, by="HGNC_symbol")
#rename columns
colnames(temp_join)[6] <- "HGNC_symbol"
colnames(temp_join2)[6] <- "Ensembl_gene_ID"
#Reorder columns
my_column_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "summary")
temp_join <- temp_join[my_column_order]
temp_join2 <- temp_join2[my_column_order]

all_matching <- rbind(temp_join, temp_join2)
all_matching <- unique(all_matching)
#Prepare output for final merge.
dge_prioritized <- data.frame(dge_gxp_0=all_matching$summary, HGNC_symbol=all_matching$HGNC_symbol)
dge_prioritized <- unique(dge_prioritized)
#Make sure to group all tissues for the gene under one entry
dge_prioritized <- dge_prioritized %>% group_by(HGNC_symbol) %>% 
  mutate(dge_gxp = paste0(dge_gxp_0, collapse = ";")) 
dge_prioritized <- dge_prioritized[-c(1)]
#Join with the master table
combined <- merge(my_master, dge_prioritized, by="HGNC_symbol", all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "dge_gxp")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")

