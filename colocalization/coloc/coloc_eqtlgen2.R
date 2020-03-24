library(coloc)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (input file and variant name).n", call.=FALSE)}

my_gwas <- args[1]
my_eqtl <- args[2]

gwas <- read.table(my_gwas, header=T, stringsAsFactors=FALSE)
eqtl <- read.table(my_eqtl, header=T, stringsAsFactors=FALSE)
#Consider only top eQTL for each SNP - need to remove duplicate IDs.
maxAbsObs <- function(x) x[which.max(abs(x))]
eqtl_sub <- aggregate(eqtl$zscore, by = list(eqtl$rsid), maxAbsObs)
eqtl_sub2 <- aggregate(eqtl$pval, by = list(eqtl$rsid), min)
#Merge the eQTL results back in.
eqtl_all <- merge(eqtl_sub, eqtl_sub2, by=c("Group.1"))
#Rename the columns
colnames(eqtl_all) <- c("rsid", "zscore", "pval")
#Intersect GWAS and eQTL SNPs
keep <- intersect(gwas$RSID, eqtl_all$rsid)
#SNPs need to be in the same order in both data frames
eqtl_keep <- subset(eqtl_all, rsid %in% keep)
gwas_keep <- subset(gwas, RSID %in% eqtl_all$rsid)
eqtl_keep <- subset(eqtl_all, rsid %in% gwas_keep$RSID)
eqtl_keep <- eqtl_keep[order(eqtl_keep$rsid),] 
gwas_keep <- gwas_keep[order(gwas_keep$RSID),] 

#Merge MAF to the eQTL data frame.
eqtl_keep <- merge(x=eqtl_keep, y=gwas_keep[,c("RSID", "maf")], by.x="rsid", by.y="RSID")

#NB number of samples slightly varies in eQTLGen dataset depending on SNP/gene analysed. 
#14115 is an average number, reported in literature.
my.res <- coloc.abf(dataset1=list(pvalues=gwas_keep$PVAL, N=103066, snp=gwas_keep$RSID, type="cc", s=0.2245),
                    dataset2=list(pvalues=eqtl_keep$pval, N=14115, snp=eqtl_keep$rsid, type="quant"),
                    MAF=gwas_keep$maf)
#We are missing sdY, beta, varbeta for our eQTL dataset - ideally, would want to include those, as well for GWAS.
#order by decreasing p-value and save results to file
total_p <- my.res$summary
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
output_file <-  paste(file_path_sans_ext(my_eqtl), ".colocp", sep="")
output_file2 <-  paste(file_path_sans_ext(my_eqtl), ".totalp", sep="")
write.table(coloc_results, file=output_file, quote=F, sep="\t", row.names=F)
write.table(total_p, file=output_file2, quote=F, sep="\t", col.names=T)

#Now, for the analysis using estimates of beta and varbeta for eQTL results (single point values not available).
#Calculate standard deviation for GWAS, using formula from doi:10.1038/ng.3538.
eqtl_keep$variance <- (1 / sqrt(2 * eqtl_keep$maf * (1-eqtl_keep$maf) * (14115 + eqtl_keep$zscore^2)))^2
eqtl_keep$beta <- eqtl_keep$zscore / sqrt(2 * eqtl_keep$maf * (1-eqtl_keep$maf) * (14115 + eqtl_keep$zscore^2))
gwas_keep$variance <- (gwas_keep$SE)^2
my.res2 <- coloc.abf(dataset1=list(beta=gwas_keep$BETA, varbeta=gwas_keep$variance,N=103066,type="cc", s=0.2245),
                     dataset2=list(beta=eqtl_keep$beta, varbeta=eqtl_keep$variance,N=14115,type="quant"),
                     MAF=gwas_keep$maf)
total_p2 <- my.res2$summary
coloc_results2 <- my.res2$results
coloc_results2 <- coloc_results2[order(coloc_results2$SNP.PP.H4,decreasing=T),]
output_file3 <-  paste(file_path_sans_ext(my_eqtl), ".colocb", sep="")
output_file4 <-  paste(file_path_sans_ext(my_eqtl), ".totalb", sep="")
write.table(coloc_results2, file=output_file3, quote=F, sep="\t", row.names=F)
write.table(total_p2, file=output_file4, quote=F, sep="\t", col.names=T)
