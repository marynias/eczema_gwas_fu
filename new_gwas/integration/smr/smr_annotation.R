library("dplyr")
library("stringr")
library("tools")
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}


my_master_file <- args[1]
output <- args[2]

#2113 genes tested
pvalue_threshold = 0.05 / 2113

#Read in all the results files for individual tissues
temp = list.files(pattern="*txt", recursive = F)
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
#Heidi pvalue > 0.05 (no heterogeneity)
final_df <- final_df[final_df$p_SMR < pvalue_threshold,]
final_df <- final_df[final_df$p_HEIDI > 0.05,]
#Prepare output for final merge.
smr_prioritized <- unique(data.frame(Ensembl_gene_ID=final_df$probeID, smr0=final_df$tissue))
#Make sure to group all tissues for the gene under one entry
smr_prioritized <- smr_prioritized %>% group_by(Ensembl_gene_ID) %>% 
  mutate(smr = paste0(smr0, collapse = ";")) 
smr_prioritized <- smr_prioritized[-c(2)]
#Join with the master table
combined <- merge(my_master, smr_prioritized, by="Ensembl_gene_ID", all.x=T)
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "smr")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")


