library("dplyr")
library("stringr")
library("tools")
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
output <- args[2]

#Read in all the results files for individual tissues
temp = list.files(pattern="*tsv", recursive = F)
#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

loadFile <- function(x) {
  all_x <- str_replace(x, "SMR_", "")
  all_x <- str_replace(all_x, ".txt", "")
  print (all_x)
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  df$tissue <- all_x
  df
}

all_normal <- lapply(temp, loadFile)
final_df = as_tibble(do.call(rbind, all_normal))
#Filter out association. eQTLs p-values are already filtered to be <5x10-8.
#SMR p-vlaue has to be lower than 5x10-8 and Heidi pvalue > 0.05 (no heterogeneity)
final_df <- final_df[final_df$p_SMR < 5e-10,]
final_df <- final_df[final_df$p_HEIDI > 0.05,]
#Prepare output for final merge.
smr_prioritized <- unique(data.frame(HGNC_symbol=final_df$Gene, smr0=final_df$tissue))
#Make sure to group all tissues for the gene under one entry
smr_prioritized <- smr_prioritized %>% group_by(HGNC_symbol) %>% 
  mutate(smr = paste0(smr0, collapse = ";")) 
smr_prioritized <- smr_prioritized[-c(2)]
#Join with the master table
combined <- merge(my_master, smr_prioritized, by="HGNC_symbol", all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "smr")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")


