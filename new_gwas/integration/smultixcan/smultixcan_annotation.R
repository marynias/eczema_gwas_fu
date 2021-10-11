library("dplyr")
library("stringr")
library("tools")
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}

my_master_file <- args[1]
output <- args[2]
#Read in all the results files for individual tissues
temp = list.files(pattern="*.txt", recursive = F)
#My annotation master table
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

loadFile <- function(x) {
  all_x <- str_replace(x, "Pheno__", "")
  all_x <- str_replace(all_x, ".results.csv.txt", "")
  print (all_x)
  df <- read.csv(x, header=T, stringsAsFactors=F,row.names=NULL)
  df$tissue <- all_x
  df$alpha <- 0.05/dim(df)[1]
  #Replace ENSG identifiers
  df$gene<- str_replace(df$gene,'\\.[0-9]+' , "")
  df
}

all_normal <- lapply(temp, loadFile)
final_df = as_tibble(do.call(rbind, all_normal))

#Filter with p-values below the Bonferroni-adjusted alpha thershold.
final_df <- final_df %>% filter(pvalue < alpha)
#Prepare output for final merge.
smultixcan_prioritized <- unique(data.frame(HGNC_symbol=final_df$gene_name, smultixcan_0=final_df$tissue, Ensembl_gene_ID=final_df$gene))
#Make sure to group all tissues for the gene under one entry
smultixcan_prioritized <- smultixcan_prioritized %>% group_by(HGNC_symbol) %>% 
  mutate(smultixcan = paste0(smultixcan_0, collapse = ";")) 
smultixcan_prioritized <- smultixcan_prioritized[-c(2)]
#Join with the master table
combined <- merge(my_master, smultixcan_prioritized, by="Ensembl_gene_ID", all.x=T)
combined <- unique(combined[-c(6)])
colnames(combined)[5] <- "HGNC_symbol"
#Remove duplicate rows. 
combined <- unique(combined)
#Resort the table.
combined <- combined[gtools::mixedorder(combined$HGNC_symbol), ]
my_col_order <- c("rsid", "cytoband", "Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "smultixcan")
combined <- combined[my_col_order]
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")

