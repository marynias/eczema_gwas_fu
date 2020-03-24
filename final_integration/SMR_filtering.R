library(dplyr)

#Read in locus names
locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t")
colnames(locus_names)[2] <- "locus"

#Read in Model 13 gene results and merge with Locus names
model13 <- read.delim("Model13_gene_ranked.txt", header=T, sep="\t")
merged <- merge(model13, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=TRUE)

processed_files <- list.files("input", recursive=F, full.names=T)
#Read in results
loadFile <- function(x) {
  print (x)
  df <- read.delim(x, header=T, stringsAsFactors=F,sep="\t")
  #Remove trailing whitespace
  df$tissue <- x
  df$tissue <- gsub("input/", "", df$tissue)
  df$tissue <- gsub(".hugo", "", df$tissue)
  df
}

filter_p_val_SMR <- function(x) {
  threshold <- 0.05/length(x$probeID)
  x <- x[x$p_SMR <= threshold,]
}

filter_p_val_SMR_HEIDI <- function(x) {
  threshold <- 0.05/length(x$probeID)
  x <- x[x$p_HEIDI >= threshold,]
}

all_normal <- lapply(processed_files, loadFile)
SMR_pval_filtered <- lapply(all_normal, filter_p_val_SMR )
HEIDI_pval_filtered <- lapply(SMR_pval_filtered, filter_p_val_SMR_HEIDI)
all_normal_together <- do.call(rbind,all_normal)
all_normal_together <- as.data.frame(all_normal_together)
SMR_together <- do.call(rbind, SMR_pval_filtered)
SMR_together <- as.data.frame(SMR_together)

HEIDI_together <- do.call(rbind, HEIDI_pval_filtered)
HEIDI_together <- as.data.frame(HEIDI_together)
#20 unique genes:
length(unique(HEIDI_together$probeID))
#29 unique genes:
length(unique(SMR_together$probeID))

#Print to file
write.table(SMR_together, "all_SMR_filtered.txt", sep="\t", col.names=T, row.names=F, quote=F)
write.table(HEIDI_together, "all_HEIDI_SMR_filtered.txt", sep="\t", col.names=T, row.names=F, quote=F)

#Annotate with Model13 gene results.
all_merged <- merge(all_normal_together, merged[c("index_SNP_rsid", "locus", "rank", "gene_name")], by.x="probeID", by.y="gene_name", all.x=T)
write.table(all_merged, "all_filtered_merged.txt", sep="\t", col.names=T, row.names=F, quote=F)
SMR_merged <- merge(SMR_together, merged[c("index_SNP_rsid", "locus", "rank", "gene_name")], by.x="probeID", by.y="gene_name", all.x=T)
write.table(SMR_merged, "all_SMR_filtered_merged.txt", sep="\t", col.names=T, row.names=F, quote=F)
HEIDI_merged <- merge(HEIDI_together, merged[c("index_SNP_rsid", "locus", "rank", "gene_name")], by.x="probeID", by.y="gene_name", all.x=T)
write.table(HEIDI_merged, "all_HEIDI_filtered_merged.txt", sep="\t", col.names=T, row.names=F, quote=F)
top <- merged[merged$rank < 4,]
top_merged <- merge(all_normal_together, top[c("index_SNP_rsid", "locus", "rank", "gene_name")], by.x="probeID", by.y="gene_name", all.y=T)
write.table(top_merged, "top_filtered_merged.txt", sep="\t", col.names=T, row.names=F, quote=F)
