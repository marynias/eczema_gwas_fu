library("tidyverse")
my_gwas_p_vals <- read.table("index_snp_pval.txt", header=TRUE)
coloc_pph4 <- read.table("pph4_3mbp_all_combined.txt", header=TRUE)

merged <- merge(coloc_pph4, my_gwas_p_vals, by="RSID")
merged$tissue <- as.factor(merged$Tissue)

merged$PPH4 <- as.numeric(merged$PPH4)
merged$PVAL <- as.numeric(merged$PVAL)
merged$logPVAL <- as.numeric(log10(merged$PVAL))

ggplot(merged, aes(x=logPVAL, y=PPH4, color=Tissue)) + geom_point() 
+ geom_smooth(method=lm, se=FALSE)
