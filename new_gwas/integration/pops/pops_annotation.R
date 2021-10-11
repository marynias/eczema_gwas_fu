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
#Read in pops result file.
my_pops <- read.delim(my_significant_file, stringsAsFactors = F, header=T)
#Temprorary join of the two files by ENSG
temp_join <- merge(my_master, my_pops, by.x="Ensembl_gene_ID", by.y="ENSG_ID")
#Add missing dummy column
temp_join$ENSG_ID <- temp_join$Ensembl_gene_ID
#Temprorary join of the two files by symbol
temp_join2 <- merge(my_master, my_pops, by.x="HGNC_symbol", by.y="Symbol")
#Add missing dummy column
temp_join2$Symbol <- temp_join2$HGNC_symbol
#Rearrange columns in temp_join2 to match order in temp_join
col_order <- colnames(temp_join)
temp_join2_rearr <- temp_join2[col_order]
#Join the two temporary join
all_matching <- rbind(temp_join, temp_join2_rearr)
#Get unique rows
all_matching <- unique(all_matching)
#Remove duplicates in case of multiple scores per gene - select highest.
all_matching_temp <- as_tibble(all_matching) %>% group_by(rsid, HGNC_symbol)%>% slice_max(order_by = Score, n = 1) %>% ungroup()
#Get top 3 scores per locus
all_matching2 <- as_tibble(all_matching_temp) %>% group_by(rsid) %>% mutate(rank = dense_rank(desc(Score))) %>% slice_max(order_by = Score, n = 3)
#Prepare output for final merge.
#Combine into a df with gene names and yes as values
pops_prioritized <- data.frame(pops_prioritization_rank=all_matching2$rank, HGNC_symbol=all_matching2$HGNC_symbol)
#Join with the master table
combined <- merge(my_master, pops_prioritized, by="HGNC_symbol", all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "pops_prioritization_rank")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")


