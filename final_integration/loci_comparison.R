library(ggplot2)
#read in gene scores
gene_stats <- read.delim("Model13_gene_ranked.txt", header=T, stringsAsFactors = FALSE)

#read in file with locus-top SNP mapping
locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t", stringsAsFactors = FALSE)
colnames(locus_names)[2] <- "locus"

#Limit ourselves to loci in the paper only.
#One file per locus (save as PNG and PDF).
gene_stats_merged <- merge(gene_stats, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=T)

paper_loci <- c("1q21.3a", "1q21.3b", "2p13.3", "2q12.1", "2q37.1", "4q27", "5p13.2", "5q31.1", "6p21.32", "6p21.33", "8q21.13", "10p15.1", "10q21.2", "11p13", "11q13.1", "11q13.5", "11q24.3", "12q15", "14q13.2", "14q32.32", "16p13.13", "17q21.2", "17q25.3", "19p13.2", "20q13.33")

#keep only the 25 loci for the paper
gene_stats_merged <- gene_stats_merged[gene_stats_merged$locus %in% paper_loci,]

my_stat_table <- lapply(unique(gene_stats_merged$index_SNP_rsid), function(x) {
  my_subset <- gene_stats_merged[gene_stats_merged$index_SNP_rsid==x,]
  my_subset <- my_subset[order(-my_subset$final_score),]
  my_subset_top_10 <- my_subset[my_subset$rank<11,]
  my_subset_top_10 <- my_subset_top_10[order(-my_subset_top_10$final_score),]
  my_total <- sum(my_subset$final_score)
  my_total_top10 <- sum(my_subset_top_10$final_score)
  
  absolute_cum_sum_top10 <- cumsum(my_subset_top_10$final_score)
  percentage_cum_sum_top10 <- absolute_cum_sum_top10/my_total_top10
  absolute_cum_sum <- cumsum(my_subset_top_10$final_score)
  percentage_cum_sum <- absolute_cum_sum/my_total
  my_labels <- seq(1,10, by=1)
  plot_df <- as.data.frame(cbind(my_labels, absolute_cum_sum, percentage_cum_sum))
  plot_df_top10 <- as.data.frame(cbind(my_labels, absolute_cum_sum_top10, percentage_cum_sum_top10))
  locus <- my_subset$locus[1]
  
  my_cum_sum_plot <-  ggplot(plot_df, aes(x=my_labels,y=percentage_cum_sum)) + geom_line(size=1, color="lightgray") + geom_point(color="blue", size=2) + theme_bw() + xlab("cumulative score sum at rank x") + ylab("fraction total score") + scale_x_discrete(limits=seq(1,10)) + scale_y_continuous(breaks=seq(0,1, by=0.1), labels=seq(0,1, by=0.1), limits=c(0,1), sec.axis = sec_axis(~.*my_total, name = "total score")) + ggtitle(locus)
  
  my_cum_sum_plot_top10 <-  ggplot(plot_df_top10, aes(x=my_labels,y=percentage_cum_sum_top10)) + geom_line(size=1, color="lightgray") + geom_point(color="blue", size=2) + theme_bw() + xlab("cumulative score sum at rank x") + ylab("fraction total top 10 score") + scale_x_discrete(limits=seq(1,10)) + scale_y_continuous(breaks=seq(0,1, by=0.1), labels=seq(0,1, by=0.1), limits=c(0,1), sec.axis = sec_axis(~.*my_total_top10, name = "total top 10 score")) + ggtitle(locus)
  
  output_file_pdf <- paste(locus, "_" , x,  "_score_comp.pdf", sep="")
  output_file_png <- paste(locus, "_" , x,  "_score_comp.png", sep="")
  output_file_pdf_top10 <- paste(locus, "_" , x,  "_score_comp_top10.pdf", sep="")
  output_file_png_top10 <- paste(locus, "_" , x,  "_score_comp_top10.png", sep="")
  ggsave(output_file_pdf, plot=my_cum_sum_plot)
  ggsave(output_file_png, plot=my_cum_sum_plot)
  ggsave(output_file_pdf_top10, plot=my_cum_sum_plot_top10)
  ggsave(output_file_png_top10, plot=my_cum_sum_plot_top10)
  
  locus_label <-  rep(locus, times = 10)
  snp_label <- rep(x, times = 10)
  my_stats <- as.data.frame(cbind(locus_label, snp_label, my_labels, absolute_cum_sum_top10, percentage_cum_sum_top10))
  my_stats
})

all_stat_together <- do.call(rbind,my_stat_table)
all_stat_together <- as.data.frame(all_stat_together)
write.table(all_stat_together, file="top10_cumulative_scores.txt", sep="\t", quote=F, row.names=F)
