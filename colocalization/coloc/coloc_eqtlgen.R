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
eqtl_sub <- aggregate(eqtl$pval, by = list(eqtl$rsid), min)
#Intersect GWAS and eQTL SNPs
keep <- intersect(gwas$rsID, eqtl_sub$Group.1)
#SNPs need to be in the same order in both data frames
eqtl_keep <- subset(eqtl_sub, Group.1 %in% keep)
eqtl_keep <- eqtl_keep[order(eqtl_keep$Group.1),] 
gwas_keep <- subset(gwas, rsID %in% eqtl$rsid)
gwas_keep <- gwas_keep[order(gwas_keep$rsID),] 

#NB number of samples slightly varies in eQTLGen dataset depending on SNP/gene analysed. 
#14115 is an average number, reported in literature.
my.res <- coloc.abf(dataset1=list(pvalues=gwas_keep$p.value, N=103066, snp=gwas_keep$rsID, type="cc", s=0.2245),
                    dataset2=list(pvalues=eqtl_keep$x, N=14115, snp=eqtl_keep$Group.1, type="quant"),
                    MAF=gwas_keep$maf)
#We are missing sdY, beta, varbeta for our eQTL dataset - ideally, would want to include those, as well for GWAS.
#order by decreasing p-value and save results to file
total_p <- my.res$summary
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
output_file <-  paste(file_path_sans_ext(my_gwas), ".coloc", sep="")
output_file2 <-  paste(file_path_sans_ext(my_gwas), ".totalp", sep="")
write.table(coloc_results, file=output_file, quote=F, sep="\t", row.names=F)
write.table(total_p, file=output_file2, quote=F, sep="\t", col.names=T)

#NB same results obtained when switching dataset 1 and dataset 2 around.
#my.res2 <- coloc.abf(dataset1=list(pvalues=eqtl_keep$x, N=14115, snp=eqtl_keep$Group.1, type="quant"),
#                    dataset2=list(pvalues=gwas_keep$p.value, N=100000, snp=gwas_keep$rsID, type="cc", s=0.183),
#                    MAF=gwas_keep$eaf)

#coloc_results2 <- my.res2$results
#coloc_results2 <- coloc_results2[order(coloc_results2$SNP.PP.H4,decreasing=T),]

#COmpare results using different N sizes for GWAS results (with/without 23andme)
#Obtain very similar results (differences at 0.0001 level in terms of p-values)
#my.res3 <- coloc.abf(dataset1=list(pvalues=gwas_keep$p.value, N=40835, snp=gwas_keep$rsID, type="cc", s=0.359037508),
#                    dataset2=list(pvalues=eqtl_keep$x, N=14115, snp=eqtl_keep$Group.1, type="quant"),
#                    MAF=gwas_keep$eaf)
#my.res3$summary
#coloc_results3 <- my.res3$results
#coloc_results3 <- coloc_results3[order(coloc_results3$SNP.PP.H4,decreasing=T),]
