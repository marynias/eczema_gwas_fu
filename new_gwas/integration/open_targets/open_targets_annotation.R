library("dplyr")
library("stringr")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
input <- args[2]
output <- args[3]


#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

#My Open Targets V2G results file
final_df <- read.delim(input, header=T, stringsAsFactors=F, row.names=NULL, sep="\t")

final_df2 <- final_df %>% group_by(rsid) %>% arrange(final_table.overallScore, .by_group = TRUE) %>% mutate(rank = dplyr::dense_rank(desc(final_table.overallScore))) %>% top_n(3, final_table.overallScore)
#Prepare output for final merge.
ot_prioritized <- data.frame(open_targets_prioritization_rank=final_df2$rank, HGNC_symbol=final_df2$symbol, rsid=final_df2$rsid)
#Join with the master table
combined <- merge(my_master, ot_prioritized, by=c("HGNC_symbol", "rsid"), all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "open_targets_prioritization_rank")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")
