library(tools)
args = commandArgs(trailingOnly=TRUE)
tissue <- args[1]
interval <- args[2]

snps_file <- paste(tissue, ".snps_vs_", interval, sep="")
genes_file <- paste(tissue, ".genes", sep="")

my_snps_comm <- read.delim(snps_file)
my_genes_comm <- read.delim(genes_file)

#Merge the tables.
my_snps_all <- merge(my_snps_comm, my_genes_comm, by="COMMUNITY")
write.table(my_snps_all, file=paste(tissue, ".snps_vs_", interval, ".community", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)
