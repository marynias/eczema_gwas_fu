library("ggplot2")
library("tidyr")
interval_stats <- read.delim("onek_r2.txt", header=T, stringsAsFactors = FALSE)
interval_stats <- interval_stats[interval_stats$interval==5000,]
interval_stats <- interval_stats[!interval_stats$snp=="rs145809981",]

locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(locus_names)[2] <- "locus"
interval_stats_merged <- merge(interval_stats, locus_names, by.x="snp", by.y="Top.SNP", all.x=T)
interval_stats <- interval_stats_merged[c("locus", "length_k_1k", "all_snps_1k")]

#Filter to keep only the loci within the GWAS
interval_stats_long <- gather(interval_stats, "Parameter", "Value", c("length_k_1k", "all_snps_1k"))

my_loci <- c("1q21.3a", "1q21.3b", "2p13.3", "2q12.1", "2q37.1", "4q27", "5p13.2", "5q31.1", "6p21.32", "6p21.33", "8q21.13", "10p15.1", "10q21.2", "11p13", "11q13.1", "11q13.5", "11q24.3", "12q15", "14q13.2", "14q32.32", "16p13.13", "17q21.2", "17q25.3", "19p13.2", "20q13.33")
interval_stats_long <- interval_stats_long[interval_stats_long$locus %in% my_loci,]

#Basic stats for figure:
summary(interval_stats_long[interval_stats_long$Parameter=="length_k_1k",])
summary(interval_stats_long[interval_stats_long$Parameter=="all_snps_1k",])

interval_stats_long[interval_stats_long$Parameter=="length_k_1k",]$Value <- interval_stats_long[interval_stats_long$Parameter=="length_k_1k",]$Value / 1000 * 11985
interval_stats_long$Parameter <- as.factor(interval_stats_long$Parameter)


basic <- ggplot(data=interval_stats_long, aes(x=interval_stats_long$locus,y=interval_stats_long$Value,colour=interval_stats_long$Parameter)) + geom_point(size=4) + scale_colour_manual(values=c("purple", "green"), labels=c("SNP number", "Interval length")) + ylab(paste("Number of SNPs"))  + scale_y_continuous(sec.axis = sec_axis(~./11985, name = "Interval length [Mbp]")) + xlab(element_blank()) + scale_x_discrete(limits=my_loci)
dotplot1 <- basic + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=14), plot.title = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.title=element_blank(), legend.background=element_blank(),legend.position="top")


#227 scaling factor for log10 number of SNPs for interval length axis to be primary, and on linear scale)