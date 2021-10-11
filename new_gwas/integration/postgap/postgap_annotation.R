library("dplyr")
library("gtools")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
my_significant_file <- args[2]
output <- args[3]



#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)
#Read in POSTGAP result file.
my_postgap <- read.delim(my_significant_file, stringsAsFactors = F, header=T)
#Keep only necessary columns
my_postgap <- my_postgap[c("gene_symbol", "gene_id", "score")]
#Get unique rows
my_postgap <- unique(my_postgap)
#Temprorary join of the two files by ENSG
temp_join <- merge(my_master, my_postgap, by.x="Ensembl_gene_ID", by.y="gene_id")
#Add missing dummy column
temp_join$gene_id <- temp_join$Ensembl_gene_ID
#Temprorary join of the two files by symbol
temp_join2 <- merge(my_master, my_postgap, by.x="HGNC_symbol", by.y="gene_symbol")
#Add missing dummy column
temp_join2$gene_symbol <- temp_join2$HGNC_symbol
#Rearrange columns in temp_join2 to match order in temp_join
col_order <- colnames(temp_join)
temp_join2_rearr <- temp_join2[col_order]
#Join the two temporary join
all_matching <- rbind(temp_join, temp_join2_rearr)
#Remove duplicates in case of multiple scores per gene - select highest.
all_matching_temp <- as_tibble(all_matching) %>% group_by(rsid, gene_symbol)%>% slice_max(order_by = score, n = 1) %>% ungroup()
#Get top 3 scores per locus
all_matching2 <- all_matching_temp %>% group_by(rsid) %>% mutate(rank = dense_rank(desc(score))) %>% slice_max(order_by = score, n = 3)
#Prepare output for final merge.
#Combine into a df with gene names and yes as values
postgap_prioritized <- data.frame(POSTGAP_prioritization_rank=all_matching2$rank, HGNC_symbol=all_matching2$HGNC_symbol)
#Join with the master table
combined <- merge(my_master, postgap_prioritized, by="HGNC_symbol", all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "POSTGAP_prioritization_rank")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")
