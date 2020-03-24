
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
eqtl <- read.table(my_eqtl, header=F, stringsAsFactors=FALSE)
colnames(eqtl) <- c("gene", "rsid", "unknown", "pvalue", "zscore")  

#Consider only top eQTL for each gene - need to remove duplicate ID due to multiple gene isoforms used in the eQTL analysis.
eqtl_sub <- aggregate(eqtl$pvalue, by = list(eqtl$rsid), min)
keep <- intersect(gwas$RSID, eqtl_sub$Group.1)
#SNPs need to be in the same order in both data frames
eqtl_keep <- subset(eqtl_sub, Group.1 %in% keep)
eqtl_keep <- eqtl_keep[order(eqtl_keep$Group.1),] 
gwas_keep <- subset(gwas, RSID %in% keep)
gwas_keep <- gwas_keep[order(gwas_keep$RSID),] 
#Merge MAF to the eQTL data frame.
colnames(eqtl_keep) <- c("rsid", "pvalue") 
eqtl_keep <- merge(x=eqtl_keep, y=gwas_keep[c("RSID","maf")], by.x="rsid", by.y="RSID")

my.res <- coloc.abf(dataset1=list(pvalues=gwas_keep$PVAL, N=103066, snp=gwas_keep$RSID, type="cc", s=0.2245),
                    dataset2=list(pvalues=eqtl_keep$pvalue, N=as.integer(no_samples), snp=eqtl_keep$rsid, type="quant"),
                    MAF=gwas_keep$maf)

total_p <- my.res$summary
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
output_file <-  paste(my_output, ".colocp", sep="")
output_file2 <-  paste(my_output, ".totalp", sep="")
write.table(coloc_results, file=output_file, quote=F, sep="\t", row.names=F)
write.table(total_p, file=output_file2, quote=F, sep="\t", col.names=T)
