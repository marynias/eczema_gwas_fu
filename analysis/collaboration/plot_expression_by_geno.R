library("ggplot2")
#Load plink raw file with genotypes (dosage-based)
genotype <- read.table("rs7101927.raw", sep=" ", stringsAsFactors=F, header=T, row.names=1)
#Load RPKM values for EMSY.
expression <- read.table("rs7101927_EMSY.rpkm", sep="\t", stringsAsFactors=F, header=T, row.names=1)
expression2 <- t(as.data.frame(expression))
#Join the two tables.
combined <- merge(genotype,expression2,by="row.names")
#As factir
combined$rs7101927_G <- as.factor(combined$rs7101927_G)

#Plot box plot and scatterplot of genotype versus expression.
fig <- ggplot(combined, aes(x=rs7101927_G, y=ENSG00000158636, fill=rs7101927_G)) + geom_boxplot() + ylab("RPKM") + theme(axis.title.x=element_blank())

#Scatterplot
basic <- ggplot(data=combined, aes(x=combined$rs7101927_G, y=combined$ENSG00000158636)) + geom_point(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("Genotype") + ylab("RPKM") 