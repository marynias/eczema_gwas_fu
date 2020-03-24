library(tools)

args = commandArgs(trailingOnly=TRUE)
gwas <- read.table(args[1], header=T, stringsAsFactors=FALSE)

gwas$OR <- exp(gwas$BETA)
write.table(gwas, file=args[2], quote=F, sep="\t", col.names=T)