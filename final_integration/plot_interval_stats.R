library("ggplot2")
library("tidyr")
interval_stats <- read.delim("onek_r2.txt", header=T)
interval_stats <- interval_stats[interval_stats$interval==5000,]
interval_stats <- interval_stats[!interval_stats$snp=="rs145809981",]

locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t")
colnames(locus_names)[2] <- "locus"
interval_stats_merged <- merge(interval_stats, locus_names, by.x="snp", by.y="Top.SNP", all.x=T)
interval_stats <- interval_stats_merged[c("locus", "length_k_1k", "all_snps_1k")]

interval_stats_long <- gather(interval_stats, "Parameter", "Value", c("length_k_1k", "all_snps_1k"))
#Transform the Log number of SNPs value so we can display everything and use secondary axis label.
interval_stats_long[interval_stats_long$Parameter=="all_snps_1k",]$Value <- 2.44 * log10(interval_stats_long[interval_stats_long$Parameter=="all_snps_1k",]$Value)
interval_stats_long[interval_stats_long$Parameter=="length_k_1k",]$Value <- log2(interval_stats_long[interval_stats_long$Parameter=="length_k_1k",]$Value)
interval_stats_long$Parameter <- as.factor(interval_stats_long$Parameter)


basic <- ggplot(data=interval_stats_long, aes(x=interval_stats_long$locus,y=interval_stats_long$Value,colour=interval_stats_long$Parameter)) + geom_point(size=4) + scale_colour_manual(values=c("purple", "green"), labels=c("SNP number", "Interval length")) + ylab(paste("Log2 interval length [kbp]"))  + scale_y_continuous(sec.axis = sec_axis(~./2.44, name = "Log10 number of SNPs"), breaks = c(0,2,4,6,8,10)) + xlab(element_blank())
dotplot1 <- basic + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=14), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.title=element_blank(), legend.background=element_blank(),legend.position="top")


#227 scaling factor for log10 number of SNPs for interval length axis to be primary, and on linear scale)