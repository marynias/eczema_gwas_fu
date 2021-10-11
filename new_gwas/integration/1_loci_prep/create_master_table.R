library("gtools")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied", call.=FALSE)}


my_cyto_file <- args[1]
my_genes_file <- args[2]
output <- args[3]

my_cyto <- read.delim(my_cyto_file, stringsAsFactors = F, header=F)
my_genes <- read.delim(my_genes_file, stringsAsFactors = F, header=F)

my_cyto$cytoband <- paste(my_cyto$V1, my_cyto$V8, sep="")
my_cyto = subset(my_cyto, select = -c(V1, V2, V3, V5, V6, V7, V8, V9))
colnames(my_cyto)[1] <- "rsid"
my_genes = subset(my_genes, select = -c(V1, V2, V3, V5, V6, V7, V8, V9, V10, V14))
colnames(my_genes) <- c("rsid", "Ensembl_gene_ID", "HGNC_symbol", "HGNC_ID")
#Filter genes with no HGNC ID
my_genes <- my_genes[!is.na(my_genes$HGNC_ID),]
#Merge the two tables
combined <- merge(my_genes, my_cyto, by="rsid")
#Reorder column order
col_order <- c("rsid", "cytoband", "Ensembl_gene_ID",
               "HGNC_ID", "HGNC_symbol")
combined <- combined[, col_order]
combined <- combined[gtools::mixedorder(combined$cytoband), ]
#Write the master table to file
write.table(combined, output, quote=F, sep=",", row.names=F)