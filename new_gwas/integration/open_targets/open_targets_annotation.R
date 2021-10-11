library("dplyr")
library("stringr")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
output <- args[2]

#Read in all the scores files for individual variants.
temp = list.files(pattern="*txt")
#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

loadFile <- function(x) {
  all_x <- str_split(x, "-")
  my_rsid <- all_x[[1]][1]
  print (my_rsid)
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  df$rsid <- my_rsid
  df
}

all_normal <- lapply(temp, loadFile)
final_df = as_tibble(do.call(rbind, all_normal))
final_df2 <- final_df %>% group_by(rsid) %>% mutate(rank = dense_rank(desc(Overall.V2G))) %>% slice_max(order_by = Overall.V2G, n = 3) %>% ungroup()
#Prepare output for final merge.
ot_prioritized <- data.frame(open_targets_prioritization_rank=final_df2$rank, HGNC_symbol=final_df2$Gene, rsid=final_df2$rsid)
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
