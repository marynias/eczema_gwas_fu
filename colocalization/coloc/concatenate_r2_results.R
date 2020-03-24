library("ggplot2")
#Load in data for table
gwas <- read.table("gwas_p.txt", header=T, stringsAsFactors=FALSE) 
ukbiobank <- read.table("ukbiobank_r2.txt", header=T, stringsAsFactors=FALSE) 
onek <- read.table("onek_r2.txt", header=T, stringsAsFactors=FALSE) 
eyeballed <- read.table("eyeballed.txt", header=T, stringsAsFactors=FALSE)  
combined_a <- merge( gwas, ukbiobank, by=c("snp", "interval"), all=TRUE)
combined <- merge(combined_a, onek, by=c("snp", "interval"), all=TRUE)
#Merge with eyeballed
combined2 <- merge(combined, eyeballed[c("snp", "length_k")], by="snp", all=TRUE)
write.table(combined2, file="all_r2.summary", quote=F, sep="\t", row.names=F)

#Load in data for plotting
gwas2 <- read.table("gwas_p2.txt", header=T, stringsAsFactors=FALSE) 
ukbiobank2 <- read.table("ukbiobank_r22.txt", header=T, stringsAsFactors=FALSE) 
onek2 <- read.table("onek_r22.txt", header=T, stringsAsFactors=FALSE) 
all_a <- rbind(gwas2, ukbiobank2)
all <- rbind(all_a, onek2)
#Bind with eyeballed - first drop unnecessary columns
all_plot <- all[c("snp", "length_k", "method", "interval")]
#Seperate for 1Mbp and 1Mbp intervals
all_plot_1Mbp <- all_plot[all_plot$interval==5000,]
all_plot_3Mbp <-  all_plot[all_plot$interval==15000,]

all_eyeballed_1Mbp <- rbind(all_plot_1Mbp, eyeballed)
all_eyeballed_3Mbp <- rbind(all_plot_3Mbp, eyeballed)
##Size interval plot
#1Mbp
basic <- ggplot(data=all_eyeballed_1Mbp, aes(x=all_eyeballed_1Mbp$snp,y=all_eyeballed_1Mbp$length_k,colour=method)) + geom_point(size=3) + scale_colour_manual(values=c("black", "red2", "green", "cyan")) + xlab("SNP id") + ylab(paste("Interval length in kbp")) + ggtitle("Max 1Mbp around index SNP") 
dotplot1 <- basic + theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=18), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("size_interval_1Mbp2.pdf", dotplot1, dpi=300, height=8, width=20)
#3Mbp
basic2 <- ggplot(data=all_eyeballed_3Mbp, aes(x=all_eyeballed_3Mbp$snp,y=all_eyeballed_3Mbp$length_k,colour=method)) + geom_point(size=3) + scale_colour_manual(values=c("black", "red2", "green", "cyan")) + xlab("SNP id") + ylab(paste("Interval length in kbp")) + ggtitle("Max 3Mbp around index SNP") 
dotplot2 <- basic2 + theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=18), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("size_interval_3Mbp2.pdf", dotplot2, dpi=300, height=8, width=20)


##Number of SNPs plot
all_1Mbp <- all[all$interval==5000,]
all_3Mbp <-  all[all$interval==15000,]

#1Mbp
basic <- ggplot(data=all_1Mbp, aes(x=all_1Mbp$snp,y=log10(all_1Mbp$all_snps),colour=method)) + geom_point(size=3) + scale_colour_manual(values=c("red2", "green", "cyan")) + xlab("SNP id") + ylab(paste("Log10 number of SNPs")) + ggtitle("Max 1Mbp around index SNP") 
dotplot1 <- basic + theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=18), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("number_snps_1Mbp2.pdf", dotplot1, dpi=300, height=8, width=20)

#3Mbp
basic2 <- ggplot(data=all_3Mbp, aes(x=all_3Mbp$snp,y=log10(all_3Mbp$all_snps),colour=method)) + geom_point(size=3) + scale_colour_manual(values=c("red2", "green", "cyan")) + xlab("SNP id") + ylab(paste("Log10 number of SNPs")) + ggtitle("Max 3Mbp around index SNP") 
dotplot2 <- basic2 + theme(axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=18), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("number_snps_3Mbp2.pdf", dotplot1, dpi=300, height=8, width=20)