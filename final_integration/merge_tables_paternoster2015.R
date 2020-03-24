library("dplyr")
library("plyr")
library("reshape2")
library("ggplot2")
locus_discovery <- read.delim("Paternoster2015_loci_discovery.txt", header=T, sep="\t", stringsAsFactors = F)
locus_discovery <- tbl_df(locus_discovery)
locus_replication <- read.delim("Paternoster2015_loci_replication.txt", header=T, sep="\t", stringsAsFactors = F)
locus_replication <- tbl_df(locus_replication)

model13_gene <- read.delim("Model13_gene_ranked.txt", sep="\t", header=T)
model13_gene_top3 <- model13_gene[model13_gene$rank < 4,]
model13_gene_top3 <- tbl_df(model13_gene_top3)
#Read in top ranked SNP for each locus - Model09.
model09_snp <- read.delim("Model09_snp_ranked.txt", sep="\t", header=T)
model09_snp_top3 <- model09_snp[model09_snp$rank < 4,]
model09_snp_top3 <- tbl_df(model09_snp_top3)

#Read in top ranked SNP for each locus - Model10.
model10_snp <- read.delim("Model10_snp_ranked.txt", sep="\t", header=T)
model10_snp_top3 <- model10_snp[model10_snp$rank < 4,]
model10_snp_top3 <- tbl_df(model10_snp_top3)

model14_snp <- read.delim("Model14_snp_ranked.txt", sep="\t", header=T)
model14_snp_top3 <- model14_snp[model14_snp$rank < 4,]
model14_snp_top3 <- tbl_df(model14_snp_top3)

get_collapsed_gene <- function(x) {
my_locus <- data.frame(Variant= character(), top3= character())
for (i in unique(x$index_SNP_rsid))
{
new_column <- x[x$index_SNP_rsid==i,]$gene_name
collapsed <- paste(new_column, collapse=",")

temp <- data.frame(Variant=i, top3=collapsed)
my_locus <- rbind(my_locus, temp)
}
tbl_df(my_locus)
}

get_collapsed_snp <- function(x) {
  my_locus <- data.frame(Variant= character(), top3= character())
  for (i in unique(x$index_SNP_rsid))
  {
    new_column <- x[x$index_SNP_rsid==i,]$current_SNP_rsid
    collapsed <- paste(new_column, collapse=",")
    
    temp <- data.frame(Variant=i, top3=collapsed)
    my_locus <- rbind(my_locus, temp)
  }
  tbl_df(my_locus)
}

model13_gene_top3_collapsed <- get_collapsed_gene(model13_gene_top3)
model09_snp_top3_collapsed <- get_collapsed_snp(model09_snp_top3)
model10_snp_top3_collapsed <- get_collapsed_snp(model10_snp_top3)
model14_snp_top3_collapsed <- get_collapsed_snp(model14_snp_top3)

merged_gene_discovery <- join(locus_discovery, model13_gene_top3_collapsed, by="Variant", type="left")

merged_gene_replication <- join(locus_replication, model13_gene_top3_collapsed, by="Variant", type="left")

merged_snp_discovery_model09 <- join(locus_discovery, model09_snp_top3_collapsed, by="Variant", type="left")

merged_snp_replication_model09 <- join(locus_replication, model09_snp_top3_collapsed, by="Variant", type="left")

merged_snp_discovery_model10 <- join(locus_discovery, model10_snp_top3_collapsed, by="Variant", type="left")

merged_snp_replication_model10 <- join(locus_replication, model10_snp_top3_collapsed, by="Variant", type="left")

merged_snp_discovery_model14 <- join(locus_discovery, model14_snp_top3_collapsed, by="Variant", type="left")

merged_snp_replication_model14 <- join(locus_replication, model14_snp_top3_collapsed, by="Variant", type="left")

write.table(merged_gene_discovery, "merged_gene_discovery.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(merged_gene_replication, "merged_gene_replication.txt", sep="\t", quote=F, col.names=T, row.names=F)

write.table(merged_snp_discovery_model09, "merged_snp_discovery_model09.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(merged_snp_replication_model09, "merged_snp_replication_model09.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(merged_snp_discovery_model10, "merged_snp_discovery_model10.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(merged_snp_replication_model10, "merged_snp_replication_model10.txt", sep="\t", quote=F, col.names=T, row.names=F)

write.table(merged_snp_discovery_model14, "merged_snp_discovery_model14.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(merged_snp_replication_model14, "merged_snp_replication_model14.txt", sep="\t", quote=F, col.names=T, row.names=F)

#Plot p-values of the top 3 SNPs - compare model 9 to model 10.
#gwas <-  read.delim("results.euro.pval.1k.dbsnp", header=T, sep="\t", stringsAsFactors = F)
#Model10 <- merge(model10_snp_top3, gwas, by.x="current_SNP_rsid", by.y="RSID")
#Model14 <- merge(model14_snp_top3, gwas, by.x="current_SNP_rsid", by.y="RSID")
#Model09 <- merge(model09_snp_top3, gwas, by.x="current_SNP_rsid", by.y="RSID")
#rm(gwas)
#write.table(Model10, "Model10_p.value.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(Model09, "Model09_p.value.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(Model14, "Model14_p.value.txt", sep="\t", quote=F, col.names=T, row.names=F)

Model10 <- read.delim("Model10_p.value.txt", sep="\t", header=T)
Model09 <- read.delim("Model09_p.value.txt", sep="\t", header=T)
Model14 <- read.delim("Model14_p.value.txt", sep="\t", header=T)

Model10_dropped <- Model10[c("current_SNP_rsid", "PVAL")]
Model10_dropped$Model <- "Model10"
summary(Model10_dropped$PVAL)
Model09_dropped <- Model09[c("current_SNP_rsid", "PVAL")]
Model09_dropped$Model <- "Model09"
Model14_dropped <- Model14[c("current_SNP_rsid", "PVAL")]
Model14_dropped$Model <- "Model14"
summary(Model14_dropped$PVAL)


all <- rbind(Model09_dropped, Model10_dropped, Model14_dropped)

ggplot(Model10, aes(x=PVAL)) + geom_histogram(bins=20, fill="darkgreen", color="black") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))  + scale_x_continuous(trans='log10')

ggplot(Model09, aes(x=PVAL)) + geom_histogram(bins=20, fill="purple", color="black") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))  + scale_x_continuous(trans='log10')

ggplot(Model14, aes(x=PVAL)) + geom_histogram(bins=20, fill="purple", color="black") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))  + scale_x_continuous(trans='log10')

ggplot(data=all, aes(all$PVAL, fill=all$Model)) + geom_histogram(col="black", position='dodge', bindwith=1) + scale_fill_manual(values=c("black", "red2", "darkgreen"), name="SNP Model")  + scale_x_continuous(trans='log10') + xlab("p-value") 