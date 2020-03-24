library(coloc)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (input file and variant name).n", call.=FALSE)}

no_samples <- args[1]
my_output <- args[2]
my_gwas <- args[3]
my_eqtl <- args[4]

gwas <- read.table(my_gwas, header=T, stringsAsFactors=FALSE)
eqtl <- read.table(my_eqtl, header=T, stringsAsFactors=FALSE)

keep <- intersect(gwas$RSID, eqtl$rsid)
#SNPs need to be in the same order in both data frames
eqtl_keep <- subset(eqtl, rsid %in% keep)
gwas_keep <- subset(gwas, RSID %in% keep)
eqtl_keep <- eqtl_keep[order(eqtl_keep$rsid),] 
gwas_keep <- gwas_keep[order(gwas_keep$RSID),] 

#Merge MAF to the eQTL data frame.
eqtl_keep <- merge(x=eqtl_keep, y=gwas_keep[c("RSID","maf")], by.x="rsid", by.y="RSID")

my.res <- coloc.abf(dataset1=list(pvalues=gwas_keep$PVAL, N=103066, snp=gwas_keep$RSID, type="cc", s=0.2245),
                    dataset2=list(pvalues=eqtl_keep$pval_nominal, N=as.integer(no_samples), snp=eqtl_keep$rsid, type="quant"),
                    MAF=gwas_keep$maf)

total_p <- my.res$summary
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
output_file <-  paste(my_output, ".colocp", sep="")
output_file2 <-  paste(my_output, ".totalp", sep="")
write.table(coloc_results, file=output_file, quote=F, sep="\t", row.names=F)
write.table(total_p, file=output_file2, quote=F, sep="\t", col.names=T)

#Now, for the analysis using estimates of beta and varbeta for eQTL results 
gwas_keep$variance <- (gwas_keep$SE)^2
eqtl_keep$variance <- (eqtl_keep$slope_se)^2
my.res2 <- coloc.abf(dataset1=list(beta=gwas_keep$BETA, varbeta=gwas_keep$variance,N=103066,type="cc", s=0.2245),
                     dataset2=list(beta=eqtl_keep$slope, varbeta=eqtl_keep$variance,N=as.integer(no_samples), type="quant"),
                     MAF=gwas_keep$maf)
total_p2 <- my.res2$summary
coloc_results2 <- my.res2$results
coloc_results2 <- coloc_results2[order(coloc_results2$SNP.PP.H4,decreasing=T),]
output_file3 <-  paste(my_output, ".colocb", sep="")
output_file4 <-  paste(my_output, ".totalb", sep="")
print (output_file3)
write.table(coloc_results2, file=output_file3, quote=F, sep="\t", row.names=F)
write.table(total_p2, file=output_file4, quote=F, sep="\t", col.names=T)