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
#Read in coloc result file for transcripts.
#coloc_tx <- read.delim("tx_coloc_all_results.txt", stringsAsFactors = F, header=T)
#Read in coloc result file for genes.
coloc_gxp <- read.delim(my_significant_file, stringsAsFactors = F, header=T)

#Filter by mininimum PPH4 = 80%.
#Maximum PPH4 for transcripts at 0.61 at rs6419573 the IL18 locus so not worth running it in the future again.
select <- coloc_gxp[coloc_gxp$PP.H4.abf > 0.95,]
select$summary <- paste(select$tissue, " (", select$study, ")", sep="")
#Prepare output for final merge.
coloc_prioritized <- unique(data.frame(HGNC_symbol=select$hugo_name, coloc0=select$summary, rsid=select$rsid))#, Ensembl_gene_ID=select$ensembl_id))
coloc_prioritized <- unique(coloc_prioritized)
#Make sure to group all tissues for the gene under one entry
coloc_prioritized2 <- as_tibble(coloc_prioritized) %>% group_by(HGNC_symbol, rsid) %>% 
  mutate(coloc = paste0(coloc0, collapse = ";")) 
#Drop not needed column
coloc_prioritized2 <- coloc_prioritized2[-c(2)]
combined <- merge(my_master, coloc_prioritized2, by=c("HGNC_symbol", "rsid"), all.x=T)
#combined2 <- merge(my_master, coloc_prioritized2, by=c("Ensembl_gene_ID", "rsid", all.x=T))
#Resort the table.
combined <- combined[gtools::mixedorder(combined$cytoband), ]
my_col_order <- c("rsid",	"cytoband",	"Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol", "coloc")
combined <- combined[my_col_order]
combined <- unique(combined)
#Write table to file
write.table(combined, output, quote=F, sep=",", row.names=F, na="")

