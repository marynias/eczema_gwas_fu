library(ggplot2)
#read in gene scores
gene_stats <- read.delim("Model13_gene_model.txt", header=T, stringsAsFactors = FALSE)
#read in variant scores
variant_stats <- read.delim("Model14_snp_model.txt", header=T, stringsAsFactors = FALSE)
#read in file with locus-top SNP mapping
locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(locus_names)[2] <- "locus"
#Limit ourselves to loci in the paper only.

#One file per locus (save as PNG and PDF).
gene_stats_merged <- merge(gene_stats, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=T)
variant_stats_merged <- merge(variant_stats, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=T)
paper_loci <- c("1q21.3a", "1q21.3b", "2p13.3", "2q12.1", "2q37.1", "4q27", "5p13.2", "5q31.1", "6p21.32", "6p21.33", "8q21.13", "10p15.1", "10q21.2", "11p13", "11q13.1", "11q13.5", "11q24.3", "12q15", "14q13.2", "14q32.32", "16p13.13", "17q21.2", "17q25.3", "19p13.2", "20q13.33")
#Loop over each index SNP
lapply(unique(gene_stats_merged$index_SNP_rsid), function(x) {
locus <- unique(gene_stats_merged[gene_stats_merged$index_SNP_rsid==x,]$Locus)[1]
my_subset <- gene_stats_merged[gene_stats_merged$index_SNP_rsid==x,]
my_plot <- ggplot(my_subset, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "white",size=0.1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black', size = 0.2)) + scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 10)) + xlab("final gene score") + ggtitle(locus) 
output_file_pdf <- paste(locus, "_" , x,  ".pdf", sep="")
output_file_png <- paste(locus, "_" , x,  ".png", sep="")
if (locus %in% paper_loci) {
  ggsave(output_file_pdf, plot=my_plot, width = 5.27, height =4.32)
  ggsave(output_file_png, plot=my_plot, width = 5.27, height =4.32)
}
})

lapply(unique(variant_stats_merged$index_SNP_rsid), function(x) {
  locus <- unique(variant_stats_merged[variant_stats_merged$index_SNP_rsid==x,]$Locus)[1]
  my_subset <- variant_stats_merged[variant_stats_merged$index_SNP_rsid==x,]
  my_plot <- ggplot(my_subset, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "white",size=0.1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black', size = 0.2)) + scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 10)) + xlab("final SNP score") + ggtitle(locus) 
  output_file_pdf <- paste(locus, "_" , x,  ".pdf", sep="")
  output_file_png <- paste(locus, "_" , x,  ".png", sep="")
  if (locus %in% paper_loci) {
    ggsave(output_file_pdf, plot=my_plot, width = 5.27, height =4.32)
    ggsave(output_file_png, plot=my_plot, width = 5.27, height =4.32)
  }
})


