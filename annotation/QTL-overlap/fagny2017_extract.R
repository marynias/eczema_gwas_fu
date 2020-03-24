library(tools)
args = commandArgs(trailingOnly=TRUE)
tissue <- args[1]
my_file <- paste(tissue, ".rds", sep="")
my_results <- readRDS(my_file, refhook = NULL)
write.table(my_results$qtls, file=paste(tissue, ".qtls", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)

my_snps_comm <- my_results$snps
my_genes_comm <- my_results$genes
my_go_comm <- my_results$go

#Merge the tables.
write.table(my_snps_comm, file=paste(tissue, ".snps", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)
write.table(my_genes_comm, file=paste(tissue, ".genes", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)
write.table(my_go_comm, file=paste(tissue, ".go", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)