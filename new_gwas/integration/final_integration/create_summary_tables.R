library("dplyr")
library("stringr")
#Merge results tables to master table.
my_master <- read.csv("paternoster2015_master.csv", stringsAsFactors = F, header=T)

#Read in results tables and merge them
temp = list.files(pattern="paternoster2015*")
temp <- temp[!temp %in% c("paternoster2015_master.csv", "paternoster2015_vep.csv")]

my_cols = colnames(my_master)

for (my_table in temp) {
  my_table_read <- read.csv(my_table, stringsAsFactors = F, header=T)
  #Check if no error reading in VEP resulsts table - if so, fix manually by appending to the master table first in R.
  my_master <- merge(my_master, my_table_read, by=my_cols)

}

#Convert all gaps to NA.
my_master <- as_tibble(my_master) %>% mutate_all(list(~na_if(.,"")))
number_of_datasources = dim(my_master)[2] - 5
my_master$total_evidence_sources <- number_of_datasources - rowSums(is.na(my_master)) 

#Counting number of multiple entries per cell (seperated with ";")
my_master$total_evidence_pieces <- my_master$total_evidence_sources + apply(my_master,1, function(x)sum(str_count(x, ";"), na.rm=T))

#Sort in descending order and save to file.
my_master <- my_master[order(-my_master$total_evidence_sources, -my_master$total_evidence_pieces),]
write.table(my_master, "final_summary_table.csv", quote=F, sep=",", row.names=F, na="")

#Seperate out by rsid and save to file
my_rsids <- unique(my_master$rsid)

saveFile <- function(x) {
my_temp_df <- my_master[my_master$rsid == x,]  
my_temp_file <- paste(x, "_summary_table.csv", sep="")
write.table(my_temp_df, my_temp_file, quote=F, sep=",", row.names=F, na="")  
}

lapply(my_rsids, saveFile)

#Create a table for graphical summary.
fig_score <- my_master %>% group_by(rsid) %>% slice_max(order_by = total_evidence_sources, n = 1) %>% ungroup()
#Now group by evidence pieces in case of ties
fig_score <- fig_score %>% group_by(rsid) %>% slice_max(order_by = total_evidence_pieces, n = 1) %>% ungroup()
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
#Replace all NA with 0
fig_score[is.na(fig_score)] <- 0
write.table(fig_score, "figure_summary_table.csv", quote=F, sep=",", row.names=F, na="")
