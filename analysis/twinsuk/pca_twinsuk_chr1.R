library(SNPRelate)

bed.fn <- "data.chr01.bed"
fam.fn <- "data.chr01.fam"
bim.fn <- "data.chr01.bim"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "data.chr01.gds")
genofile <- snpgdsOpen("data.chr01.gds")

#LD-based SNP pruning.
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

# Get all selected snp id
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                  EV5 = pca$eigenvect[,5],    # the fifth eigenvector
                  EV6 = pca$eigenvect[,6],    # the sixth eigenvector
                  stringsAsFactors = FALSE)

write.table(tab, file="data.chr01.pca", sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

#ibd <- snpgdsIBDMoM(genofile, sample.id=pca$sample.id, snp.id=snpset.id,
#                    maf=0.05, missing.rate=0.05, num.thread=2)