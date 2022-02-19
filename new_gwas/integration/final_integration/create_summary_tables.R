library("dplyr")
library("stringr")
library("tools")
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}

my_master_file <- args[1]
lead_SNP_tab <- args[2]
gwas_name <- args[3]

#Merge results tables to master table.
my_master <- read.csv(my_master_file, stringsAsFactors = F, header=T)

#Read in results from lead SNP table. Use to merge to give GRCh37 variant locations.
lead_SNP_table <- read.delim(lead_SNP_tab, stringsAsFactors = F, header=T, sep=",")

my_master <- merge(my_master, lead_SNP_table[c("RSID", "CHR", "POS")], by.x="rsid", by.y="RSID")

my_master <- my_master %>% dplyr::select(rsid, cytoband, CHR, POS, Ensembl_gene_ID, HGNC_ID, HGNC_symbol)

#Read in results tables and merge them
my_pattern <- paste0(gwas_name, "*.csv")
temp = Sys.glob(my_pattern)

my_cols = c("rsid", "cytoband", "Ensembl_gene_ID", "HGNC_ID", "HGNC_symbol")

for (my_table in temp) {
  my_table_read <- read.csv(my_table, stringsAsFactors = F, header=T)
  #Check if no error reading in VEP resulsts table - if so, fix manually by appending to the master table first in R.
  my_master <- merge(my_master, my_table_read, by=my_cols)

}

#Convert all gaps to NA.
my_master <- as_tibble(my_master) %>% mutate_all(list(~na_if(.,"")))

number_of_datatypes = dim(my_master)[2] - 7
my_master$total_evidence_types <- number_of_datatypes - rowSums(is.na(my_master)) 

#Total evidence pieces without xQTL colocalisation-like evidence.
my_master_subset <- my_master %>% select(-smr, -coloc, -smultixcan)
#temprorary total evidence pieces (without coloc methods)
my_master$total_evidence_pieces <- number_of_datatypes - 3 - rowSums(is.na(my_master_subset)) 
my_master$total_evidence_pieces <-  my_master$total_evidence_pieces + apply(my_master_subset,1, function(x)sum(str_count(x, ";"), na.rm=T))

#Counting number of multiple entries per cell (seperated with ";")
my_master$coloc_list <- strsplit(my_master$coloc, ";")
my_master$smr_list <- strsplit(my_master$smr, ";")
my_master$smultixcan_list <- strsplit(my_master$smultixcan, ";")

for (a in 1:nrow(my_master)) {
  listoflists <- c(my_master$coloc_list[a], my_master$smr_list[a], my_master$smultixcan_list[a])
  temp <- do.call(c, listoflists)
  my_master$total_evidence_pieces[a] <- my_master$total_evidence_pieces[a] + length(unique(na.omit(temp)))
}
#Drop unncessary columns
my_master <- my_master %>% select(-coloc_list, -smr_list, -smultixcan_list)
#Sort in descending order and save to file.
my_master <- my_master[order(-my_master$total_evidence_pieces, -my_master$total_evidence_types),]
output_file <- paste("final_summary_table_", gwas_name, ".csv", sep="")
write.table(my_master, output_file, quote=F, sep=",", row.names=F, na="")

#Seperate out by rsid and save to file
my_rsids <- unique(my_master$rsid)

saveFile <- function(x) {
my_temp_df <- my_master[my_master$rsid == x,]  
my_temp_file <- paste(x, "_", gwas_name, "_summary_table.csv", sep="")
write.table(my_temp_df, my_temp_file, quote=F, sep=",", row.names=F, na="")  
}

lapply(my_rsids, saveFile)

#Create a table for graphical summary. Min. 3 evidence types
fig_score <- my_master %>% filter(total_evidence_types > 2) %>% group_by(rsid) %>% slice_max(order_by = total_evidence_pieces, n = 1) %>% ungroup()
#Now group by evidence pieces in case of ties
fig_score <- fig_score %>% group_by(rsid) %>% slice_max(order_by = total_evidence_types, n = 1) %>% ungroup()
#Prepare a table for output, which provides a count of entries in the coloc, SMR, Smultixcan, DGE proteome, DGE expression columns
fig_score$coloc <- str_count(fig_score$coloc, "\\(")
fig_score$smr <- ifelse(is.na(fig_score$smr), 0, 1+str_count(fig_score$smr, ";"))
fig_score$smultixcan <- ifelse(is.na(fig_score$smultixcan), 0, 1+str_count(fig_score$smultixcan, ";"))
fig_score$dge_gxp <- str_count(fig_score$dge_gxp, "\\(")
fig_score$dge_proteome <- str_count(fig_score$dge_proteome, "\\(")
#Change MendelVar output to binary
fig_score$MendelVar_sig_enrichment <- ifelse(fig_score$MendelVar_sig_enrichment  == "yes", 1, 0)
fig_score$MendelVar_skin_keywords <- ifelse(fig_score$MendelVar_skin_keywords  == "yes", 1, 0)
#Change DEPICT output to binary
fig_score$DEPICT_prioritization <- ifelse(fig_score$DEPICT_prioritization  == "yes", 1, 0)
#Change MAGMA output to binary
fig_score$MAGMA_prioritization <- ifelse(fig_score$MAGMA_prioritization  == "yes", 1, 0)
fig_score$VEP_intron <- ifelse(fig_score$VEP_intron  == "yes", 1, 0)
fig_score$VEP_missense <- ifelse(fig_score$VEP_missense  == "yes", 1, 0)
#Replace all NA with 0
fig_score[is.na(fig_score)] <- 0
output_file2 <- paste("figure_summary_table_", gwas_name, ".csv", sep="")
write.table(fig_score, output_file2, quote=F, sep=",", row.names=F, na="")
