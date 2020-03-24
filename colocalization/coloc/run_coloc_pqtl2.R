library(coloc)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (input file and variant name).n", call.=FALSE)}

my_output <- args[1]
my_gwas <- args[2]
my_eqtl <- args[3]

gwas <- read.table(my_gwas, header=T, stringsAsFactors=FALSE)
eqtl <- read.table(my_eqtl, header=T, stringsAsFactors=FALSE)
#Intersect GWAS and eQTL SNPs
keep <- intersect(gwas$RSID, eqtl$snp)
#SNPs need to be in the same order in both data frames
eqtl_keep <- subset(eqtl, snp %in% keep)
eqtl_keep <- eqtl_keep[order(eqtl_keep$snp),] 
gwas_keep <- subset(gwas, RSID %in% eqtl_keep$snp)
gwas_keep <- gwas_keep[order(gwas_keep$RSID),] 

#P-value based analysis.
my.res <- coloc.abf(dataset1=list(pvalues=gwas_keep$PVAL, N=103066, snp=gwas_keep$RSID, type="cc", s=0.2245),
                    dataset2=list(pvalues=eqtl_keep$p, N=3301, snp=eqtl_keep$snp, type="quant"),
                    MAF=gwas_keep$maf)
total_p <- my.res$summary
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
output_file <-  paste(my_output, ".colocp", sep="")
output_file2 <-  paste(my_output, ".totalp", sep="")
write.table(coloc_results, file=output_file, quote=F, sep="\t", row.names=F)
write.table(total_p, file=output_file2, quote=F, sep="\t", col.names=T)

#Carry out the analysis using beta and its variance.
gwas_keep$variance <- (gwas_keep$SE)^2
eqtl_keep$variance <- (eqtl_keep$se)^2
my.res2 <- coloc.abf(dataset1=list(beta=gwas_keep$BETA, varbeta=gwas_keep$variance,N=103066,type="cc", s=0.2245),
                     dataset2=list(beta=eqtl_keep$beta, varbeta=eqtl_keep$variance,N=3301,type="quant"),
                     MAF=gwas_keep$maf)
total_p2 <- my.res2$summary
coloc_results2 <- my.res2$results
coloc_results2 <- coloc_results2[order(coloc_results2$SNP.PP.H4,decreasing=T),]
output_file3 <-  paste(my_output, ".colocb", sep="")
output_file4 <-  paste(my_output, ".totalb", sep="")
write.table(coloc_results2, file=output_file3, quote=F, sep="\t", row.names=F)
write.table(total_p2, file=output_file4, quote=F, sep="\t", col.names=T)
