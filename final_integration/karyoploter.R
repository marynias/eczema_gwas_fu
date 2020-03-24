library("karyoploteR")
#Read in the coordinates for all the loci.
loci_r2 <- read.delim("interval_r2_0.2_1k_custom.bed", header=F, stringsAsFactors = F)
colnames(loci_r2) <- c("chrom", "start","end")
my_index_snps <-  read.delim("paternoster_2015_lead_SNPs.txt", header=T, stringsAsFactors = F)

#Read in the annotation.
my_annotations <- read.delim("paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches_ref.hugo", header=T, stringsAsFactors = F)

#Read in the scores.
scores <- read.delim("Model13_gene_ranked.txt", header=T, stringsAsFactors = F)

for (i in seq_along(loci_r2$chrom)) {
custom.genome <- toGRanges(loci_r2[i,])
custom_annotation <- my_annotations[my_annotations$snp_name==loci_r2[i,1],]

seq_length <- abs(loci_r2[i,3] - loci_r2[i,2])

#Subset to scores of interest
custom_scores <- scores[scores$index_SNP_rsid==loci_r2[i,1],]
#Adjust final score to top value.
max_value_score <- max(custom_scores$final_score)
custom_scores$adjusted_score <- custom_scores$final_score/max_value_score * 100
#Join with the annotation object
my_merged_annotation <- merge(custom_annotation, custom_scores, by.x="ensembl", by.y="gene_name", all.x=T)
my_merged_annotation <- my_merged_annotation[c("snp_name", "starts", "ends", "chr", "ensembl", "gene_name", "final_score", "adjusted_score", "rank")]
my_merged_annotation <- my_merged_annotation[order(my_merged_annotation$rank),]
#Drop genes with no score.
my_merged_annotation <- subset(my_merged_annotation, !is.na(final_score))
#Get middle of the gene.
my_merged_annotation$middle_gene <- round((my_merged_annotation$starts + my_merged_annotation$ends) / 2)

#Produce output for plotting with gassocplot
##Unadjusted score
my_out <- paste(loci_r2[i,1], "_genes.unadjusted.gassocplot", sep="")
my_df <- my_merged_annotation[c("gene_name", "chr", "middle_gene", "final_score")]
write.table(my_df, my_out, col.names=T, quote=F, row.names=F, sep="\t")
##Adjusted score
my_out2 <- paste(loci_r2[i,1], "_genes.adjusted.gassocplot", sep="")
my_df2 <- my_merged_annotation[c("gene_name", "chr", "middle_gene", "adjusted_score")]
write.table(my_df2, my_out2, col.names=T, quote=F, row.names=F, sep="\t")
}

##Produce output for Gassocplot for plotting of the SNPs.

#Read in annotations.
onek_annotation <- read.delim("Model09_snp_ranked.1k", header=F, sep="\t")
onek_annotation <- onek_annotation[c("V12", "V2", "V3")]
colnames(onek_annotation) <- c("rsid", "chrom_1k", "pos_1k")
dbsnp_annotation <- read.delim("Model09_snp_ranked.dbsnp", header=F, sep=" ")
dbsnp_annotation <- dbsnp_annotation[c("V1", "V6", "V7")]
colnames(dbsnp_annotation) <- c("rsid", "chrom_dbsnp", "pos_dbsnp" )
#Read in rankings
model09_snp <- read.delim("Model09_snp_ranked.txt", header=T, sep="\t")
model10_snp <- read.delim("Model10_snp_ranked.txt", header=T, sep="\t")
model14_snp <- read.delim("Model14_snp_ranked.txt", header=T, sep="\t")

#Merge
model09_snp_all <- merge(model09_snp,dbsnp_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model09_snp_all <- merge(model09_snp_all,onek_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model09_snp_all <- tbl_df(model09_snp_all)
model09_snp_all <- model09_snp_all %>% mutate(chrom=ifelse(!is.na(chrom_dbsnp), chrom_dbsnp,chrom_1k))
model09_snp_all <- model09_snp_all %>% mutate(pos=ifelse(!is.na(pos_dbsnp), pos_dbsnp,pos_1k))
#Drop SNPs with empty positions - none of them have rank higher than 300.
model09_snp_all  <- model09_snp_all %>% filter(!is.na(pos))

model10_snp_all <- merge(model10_snp,dbsnp_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model10_snp_all <- merge(model10_snp_all,onek_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model10_snp_all <- tbl_df(model10_snp_all)
model10_snp_all <- model10_snp_all %>% mutate(chrom=ifelse(!is.na(chrom_dbsnp), chrom_dbsnp,chrom_1k))
model10_snp_all <- model10_snp_all %>% mutate(pos=ifelse(!is.na(pos_dbsnp), pos_dbsnp,pos_1k))
model10_snp_all <- model10_snp_all %>%  filter(!is.na(pos))

model14_snp_all <- merge(model14_snp,dbsnp_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model14_snp_all <- merge(model14_snp_all,onek_annotation, by.x="current_SNP_rsid", by.y="rsid",all.x=T)
model14_snp_all <- tbl_df(model14_snp_all)
model14_snp_all <- model14_snp_all %>% mutate(chrom=ifelse(!is.na(chrom_dbsnp), chrom_dbsnp,chrom_1k))
model14_snp_all <- model14_snp_all %>% mutate(pos=ifelse(!is.na(pos_dbsnp), pos_dbsnp,pos_1k))
#Drop SNPs with empty positions - none of them have rank higher than 300.
model14_snp_all  <- model14_snp_all %>% filter(!is.na(pos))

#Print to file.
for (locus in unique(model10_snp_all$index_SNP_rsid))
{
temp <- model09_snp_all[model09_snp_all$index_SNP_rsid==locus,] 
temp <- temp[order(temp$pos),]
#Get adjusted score
temp$adjusted_score <- (temp$final_score / max(temp$final_score)) * 100
temp_adjusted <- temp[c("current_SNP_rsid", "chrom", "pos", "adjusted_score")]
temp_unadjusted <- temp[c("current_SNP_rsid", "chrom", "pos", "final_score")]
my_out <- paste(locus, "_snps.adjusted_Model09.gassocplot", sep="")
my_out_un <- paste(locus, "_snps.unadjusted_Model09.gassocplot", sep="")
write.table(temp_adjusted, my_out, col.names=T, quote=F, row.names=F, sep="\t")
write.table(temp_unadjusted, my_out_un, col.names=T, quote=F, row.names=F, sep="\t")

temp2 <- model10_snp_all[model10_snp_all$index_SNP_rsid==locus,]  
temp2 <- temp2[order(temp2$pos),]
temp2$adjusted_score <- (temp2$final_score / max(temp2$final_score)) * 100
temp_adjusted2 <- temp2[c("current_SNP_rsid", "chrom", "pos", "adjusted_score")]
temp_unadjusted2 <- temp2[c("current_SNP_rsid", "chrom", "pos", "final_score")]
my_out2 <- paste(locus, "_snps.adjusted_Model10.gassocplot", sep="")
my_out_un2 <- paste(locus, "_snps.unadjusted_Model10.gassocplot", sep="")
write.table(temp_adjusted2, my_out2, col.names=T, quote=F, row.names=F, sep="\t")
write.table(temp_unadjusted2, my_out_un2, col.names=T, quote=F, row.names=F, sep="\t")
}

for (locus in unique(model14_snp_all$index_SNP_rsid))
{
  temp <- model14_snp_all[model14_snp_all$index_SNP_rsid==locus,] 
  temp <- temp[order(temp$pos),]
  #Get adjusted score
  temp$adjusted_score <- (temp$final_score / max(temp$final_score)) * 100
  temp_adjusted <- temp[c("current_SNP_rsid", "chrom", "pos", "adjusted_score")]
  temp_unadjusted <- temp[c("current_SNP_rsid", "chrom", "pos", "final_score")]
  my_out <- paste(locus, "_snps.adjusted_Model14.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model14.gassocplot", sep="")
  write.table(temp_adjusted, my_out, col.names=T, quote=F, row.names=F, sep="\t")
  write.table(temp_unadjusted, my_out_un, col.names=T, quote=F, row.names=F, sep="\t")
}
  
library("gassocplot")
library("ggplot2")
library("grid")
library("gridExtra")
library("gtable")

#Remove rs145809981 synonym from the analysis
loci_r2 <- loci_r2[!loci_r2$chrom=="rs145809981",]

#Gene
for (locus in unique(loci_r2$chrom))
{
  my_out <- paste(locus, "_genes.adjusted.gassocplot", sep="")
  my_out_un <- paste(locus, "_genes.unadjusted.gassocplot", sep="")
  #Read in data.
  my_gene_df <- read.delim(my_out, sep="\t", header=T)
  my_gene_df_un <- read.delim(my_out_un, sep="\t", header=T)
  colnames(my_gene_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_gene_df_un) <- c("marker", "chr", "pos", "prob")
  plot <- assoc_plot2(my_gene_df,type="prob", ylab="score")
  plot_un <- assoc_plot2(my_gene_df_un,type="prob", ylab="score")
  plot_name <- paste(locus, "_genes.adjusted.gassocplot.png", sep="")
  plot_name_un <- paste(locus, "_genes.unadjusted.gassocplot.png", sep="")
  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)
}

#SNP
for (locus in unique(loci_r2$chrom))
{
  my_out <- paste(locus, "_snps.adjusted_Model10.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot", sep="")
  my_out2 <- paste(locus, "_snps.adjusted_Model09.gassocplot", sep="")
  my_out_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot", sep="")
  my_snp_df <- read.delim(my_out, sep="\t", header=T)
  my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)
  my_snp_df2 <- read.delim(my_out2, sep="\t", header=T)
  my_snp_df_un2 <- read.delim(my_out_un2, sep="\t", header=T)
  colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df2) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un2) <- c("marker", "chr", "pos", "prob")
  #Limit ourselves only to the interval.
  start_position <- loci_r2[loci_r2$chrom==locus,]$start
  end_position <- loci_r2[loci_r2$chrom==locus,]$end
  my_snp_df <- my_snp_df[my_snp_df$pos <= end_position & my_snp_df$pos >= start_position,]
  my_snp_df_un <- my_snp_df_un[my_snp_df_un$pos <= end_position & my_snp_df_un$pos >= start_position,]
  my_snp_df2 <- my_snp_df2[my_snp_df2$pos <= end_position & my_snp_df2$pos >= start_position,]
  my_snp_df_un2 <- my_snp_df_un2[my_snp_df_un2$pos <= end_position & my_snp_df_un2$pos >= start_position,]
  #Write filtered results to file.
  output_filename <- paste(locus, "_snps.adjusted_Model10.interval.gassocplot", sep="")
  write.table(my_snp_df, output_filename, sep="\t", col.names=T, row.names=F, quote=F)
  output_filename_un <- paste(locus, "_snps.unadjusted_Model10.interval.gassocplot", sep="")
  write.table(my_snp_df_un, output_filename_un, sep="\t", col.names=T, row.names=F, quote=F)
  
  output_filename2 <- paste(locus, "_snps.adjusted_Model09.interval.gassocplot", sep="")
  write.table(my_snp_df2, output_filename2, sep="\t", col.names=T, row.names=F, quote=F)
  output_filename_un2 <- paste(locus, "_snps.unadjusted_Model09.interval.gassocplot", sep="")
  write.table(my_snp_df_un2, output_filename_un2, sep="\t", col.names=T, row.names=F, quote=F)
  
  plot <- assoc_plot3(my_snp_df,type="prob", ylab="score")
  plot_un <- assoc_plot3(my_snp_df_un,type="prob", ylab="score")
  plot2 <- assoc_plot3(my_snp_df2,type="prob", ylab="score")
  plot_un2 <- assoc_plot3(my_snp_df_un2,type="prob", ylab="score")
  plot_name <- paste(locus, "_snps.adjusted_Model10.gassocplot.png", sep="")
  plot_name_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot.png", sep="")
  plot_name2 <- paste(locus, "_snps.adjusted_Model09.gassocplot.png", sep="")
  plot_name_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot.png", sep="")
  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)
  assoc_plot_save(plot2, plot_name2)
  assoc_plot_save(plot_un2, plot_name_un2)
  #Create a fake correlation matrix.
  #correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
  #colnames(correlation_top_snp) <- my_snp_df$marker
  #rownames(correlation_top_snp) <- my_snp_df$marker
}

#SNP with r2.
for (locus in unique(loci_r2$chrom))
  #for (locus in c("rs7512552"))
{
  my_out <- paste(locus, "_snps.adjusted_Model10.interval.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model10.interval.gassocplot", sep="")
  my_out2 <- paste(locus, "_snps.adjusted_Model09.interval.gassocplot", sep="")
  my_out_un2 <- paste(locus, "_snps.unadjusted_Model09.interval.gassocplot", sep="")
  my_snp_df <- read.delim(my_out, sep="\t", header=T)
  my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)
  my_snp_df2 <- read.delim(my_out2, sep="\t", header=T)
  my_snp_df_un2 <- read.delim(my_out_un2, sep="\t", header=T)
  colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df2) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un2) <- c("marker", "chr", "pos", "prob")
  #Read in the correlation matrix.
  my_corr <- paste(locus, "_gassocplot.ld", sep="")
  my_corr_df <- read.delim(my_corr, sep="\t", header=T, row.names=1)
  #colnames(my_corr_df) <- colnames(my_corr_df)[2:length(colnames(my_corr_df))]
  my_corr_df <- as.matrix(my_corr_df)
  plot <- assoc_plot4(my_snp_df,my_corr_df, type="prob", ylab="score")
  plot_un <- assoc_plot4(my_snp_df_un,my_corr_df, type="prob", ylab="score")
  plot2 <- assoc_plot4(my_snp_df2,my_corr_df,type="prob", ylab="score")
  plot_un2 <- assoc_plot4(my_snp_df_un2,my_corr_df,type="prob", ylab="score")
  plot_name <- paste(locus, "_snps.adjusted_Model10.gassocplot.ld.png", sep="")
  plot_name_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot.ld.png", sep="")
  plot_name2 <- paste(locus, "_snps.adjusted_Model09.gassocplot.ld.png", sep="")
  plot_name_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot.ld.png", sep="")
  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)
  assoc_plot_save(plot2, plot_name2)
  assoc_plot_save(plot_un2, plot_name_un2)
  #Create a fake correlation matrix.
  #correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
  #colnames(correlation_top_snp) <- my_snp_df$marker
  #rownames(correlation_top_snp) <- my_snp_df$marker
}



#Gene
#for (locus in c("rs7512552"))
for (locus in unique(loci_r2$chrom))
{
locus_pos = my_index_snps[my_index_snps$RSID==locus,]$POS
my_out <- paste(locus, "_genes.adjusted.gassocplot", sep="")
my_out_un <- paste(locus, "_genes.unadjusted.gassocplot", sep="")
#Read in data.
my_gene_df <- read.delim(my_out, sep="\t", header=T)
my_gene_df_un <- read.delim(my_out_un, sep="\t", header=T)
colnames(my_gene_df) <- c("marker", "chr", "pos", "prob")
colnames(my_gene_df_un) <- c("marker", "chr", "pos", "prob")
plot <- assoc_plot2a(my_gene_df, my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot_un <- assoc_plot2a(my_gene_df_un,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot_name <- paste(locus, "_genes.adjusted.gassocplot.index_snp.png", sep="")
plot_name_un <- paste(locus, "_genes.unadjusted.gassocplot_index_snp.png", sep="")
assoc_plot_save(plot, plot_name)
assoc_plot_save(plot_un, plot_name_un)
}

#SNP
for (locus in unique(loci_r2$chrom))
{
locus_pos = my_index_snps[my_index_snps$RSID==locus,]$POS
my_out <- paste(locus, "_snps.adjusted_Model10.gassocplot", sep="")
my_out_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot", sep="")
my_out2 <- paste(locus, "_snps.adjusted_Model09.gassocplot", sep="")
my_out_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot", sep="")
my_snp_df <- read.delim(my_out, sep="\t", header=T)
my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)
my_snp_df2 <- read.delim(my_out2, sep="\t", header=T)
my_snp_df_un2 <- read.delim(my_out_un2, sep="\t", header=T)
colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")
colnames(my_snp_df2) <- c("marker", "chr", "pos", "prob")
colnames(my_snp_df_un2) <- c("marker", "chr", "pos", "prob")
#Limit ourselves only to the interval.
start_position <- loci_r2[loci_r2$chrom==locus,]$start
end_position <- loci_r2[loci_r2$chrom==locus,]$end
my_snp_df <- my_snp_df[my_snp_df$pos <= end_position & my_snp_df$pos >= start_position,]
my_snp_df_un <- my_snp_df_un[my_snp_df_un$pos <= end_position & my_snp_df_un$pos >= start_position,]
my_snp_df2 <- my_snp_df2[my_snp_df2$pos <= end_position & my_snp_df2$pos >= start_position,]
my_snp_df_un2 <- my_snp_df_un2[my_snp_df_un2$pos <= end_position & my_snp_df_un2$pos >= start_position,]
#Write filtered results to file.
output_filename <- paste(locus, "_snps.adjusted_Model10.interval.gassocplot", sep="")
write.table(my_snp_df, output_filename, sep="\t", col.names=T, row.names=F, quote=F)
output_filename_un <- paste(locus, "_snps.unadjusted_Model10.interval.gassocplot", sep="")
write.table(my_snp_df_un, output_filename_un, sep="\t", col.names=T, row.names=F, quote=F)

output_filename2 <- paste(locus, "_snps.adjusted_Model09.interval.gassocplot", sep="")
write.table(my_snp_df2, output_filename2, sep="\t", col.names=T, row.names=F, quote=F)
output_filename_un2 <- paste(locus, "_snps.unadjusted_Model09.interval.gassocplot", sep="")
write.table(my_snp_df_un2, output_filename_un2, sep="\t", col.names=T, row.names=F, quote=F)

plot <- assoc_plot3a(my_snp_df,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot_un <- assoc_plot3a(my_snp_df_un,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot2 <- assoc_plot3a(my_snp_df2,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot_un2 <- assoc_plot3a(my_snp_df_un2,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
plot_name <- paste(locus, "_snps.adjusted_Model10.gassocplot.index_snp.png", sep="")
plot_name_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot.index_snp.png", sep="")
plot_name2 <- paste(locus, "_snps.adjusted_Model09.gassocplot.index_snp.png", sep="")
plot_name_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot.index_snp.png", sep="")
assoc_plot_save(plot, plot_name)
assoc_plot_save(plot_un, plot_name_un)
assoc_plot_save(plot2, plot_name2)
assoc_plot_save(plot_un2, plot_name_un2)
#Create a fake correlation matrix.
#correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
#colnames(correlation_top_snp) <- my_snp_df$marker
#rownames(correlation_top_snp) <- my_snp_df$marker
}

#SNP with r2.
for (locus in unique(loci_r2$chrom))
#for (locus in c("rs7512552"))
{
  locus_pos = my_index_snps[my_index_snps$RSID==locus,]$POS
  my_out <- paste(locus, "_snps.adjusted_Model10.interval.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model10.interval.gassocplot", sep="")
  my_out2 <- paste(locus, "_snps.adjusted_Model09.interval.gassocplot", sep="")
  my_out_un2 <- paste(locus, "_snps.unadjusted_Model09.interval.gassocplot", sep="")
  my_snp_df <- read.delim(my_out, sep="\t", header=T)
  my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)
  my_snp_df2 <- read.delim(my_out2, sep="\t", header=T)
  my_snp_df_un2 <- read.delim(my_out_un2, sep="\t", header=T)
  colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df2) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un2) <- c("marker", "chr", "pos", "prob")
  #Read in the correlation matrix.
  my_corr <- paste(locus, "_gassocplot.ld", sep="")
  my_corr_df <- read.delim(my_corr, sep="\t", header=T, row.names=1)
  #colnames(my_corr_df) <- colnames(my_corr_df)[2:length(colnames(my_corr_df))]
  my_corr_df <- as.matrix(my_corr_df)
  plot <- assoc_plot4a(my_snp_df,my_corr_df, type="prob",my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  plot_un <- assoc_plot4a(my_snp_df_un,my_corr_df, type="prob", my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  plot2 <- assoc_plot4a(my_snp_df2,my_corr_df,type="prob", my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  plot_un2 <- assoc_plot4a(my_snp_df_un2,my_corr_df,type="prob", my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  plot_name <- paste(locus, "_snps.adjusted_Model10.gassocplot.ld.index_snp.png", sep="")
  plot_name_un <- paste(locus, "_snps.unadjusted_Model10.gassocplot.ld.index_snp.png", sep="")
  plot_name2 <- paste(locus, "_snps.adjusted_Model09.gassocplot.ld.index_snp.png", sep="")
  plot_name_un2 <- paste(locus, "_snps.unadjusted_Model09.gassocplot.ld.index_snp.png", sep="")
  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)
  assoc_plot_save(plot2, plot_name2)
  assoc_plot_save(plot_un2, plot_name_un2)
  #Create a fake correlation matrix.
  #correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
  #colnames(correlation_top_snp) <- my_snp_df$marker
  #rownames(correlation_top_snp) <- my_snp_df$marker
}

for (locus in unique(loci_r2$chrom))
{
  locus_pos = my_index_snps[my_index_snps$RSID==locus,]$POS
  my_out <- paste(locus, "_snps.adjusted_Model14.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model14.gassocplot", sep="")
  my_snp_df <- read.delim(my_out, sep="\t", header=T)
  my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)
  colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")
  #Limit ourselves only to the interval.
  start_position <- loci_r2[loci_r2$chrom==locus,]$start
  end_position <- loci_r2[loci_r2$chrom==locus,]$end
  my_snp_df <- my_snp_df[my_snp_df$pos <= end_position & my_snp_df$pos >= start_position,]
  my_snp_df_un <- my_snp_df_un[my_snp_df_un$pos <= end_position & my_snp_df_un$pos >= start_position,]
  #Write filtered results to file.
  output_filename <- paste(locus, "_snps.adjusted_Model14.interval.gassocplot", sep="")
  write.table(my_snp_df, output_filename, sep="\t", col.names=T, row.names=F, quote=F)
  output_filename_un <- paste(locus, "_snps.unadjusted_Model14.interval.gassocplot", sep="")
  write.table(my_snp_df_un, output_filename_un, sep="\t", col.names=T, row.names=F, quote=F)
  
  plot <- assoc_plot3a(my_snp_df,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
  plot_un <- assoc_plot3a(my_snp_df_un,my_locus_lab=locus, my_locus_pos=locus_pos, type="prob", ylab="score")
  
  plot_name <- paste(locus, "_snps.adjusted_Model14.gassocplot.index_snp.png", sep="")
  plot_name_un <- paste(locus, "_snps.unadjusted_Model14.gassocplot.index_snp.png", sep="")

  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)
  #Create a fake correlation matrix.
  #correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
  #colnames(correlation_top_snp) <- my_snp_df$marker
  #rownames(correlation_top_snp) <- my_snp_df$marker
}

#SNP with r2.
for (locus in unique(loci_r2$chrom))
  #for (locus in c("rs7512552"))
{
  locus_pos = my_index_snps[my_index_snps$RSID==locus,]$POS
  my_out <- paste(locus, "_snps.adjusted_Model14.interval.gassocplot", sep="")
  my_out_un <- paste(locus, "_snps.unadjusted_Model14.interval.gassocplot", sep="")

  my_snp_df <- read.delim(my_out, sep="\t", header=T)
  my_snp_df_un <- read.delim(my_out_un, sep="\t", header=T)

  colnames(my_snp_df) <- c("marker", "chr", "pos", "prob")
  colnames(my_snp_df_un) <- c("marker", "chr", "pos", "prob")

  #Read in the correlation matrix.
  my_corr <- paste(locus, "_gassocplot.ld", sep="")
  my_corr_df <- read.delim(my_corr, sep="\t", header=T, row.names=1)
  #colnames(my_corr_df) <- colnames(my_corr_df)[2:length(colnames(my_corr_df))]
  my_corr_df <- as.matrix(my_corr_df)
  plot <- assoc_plot4a(my_snp_df,my_corr_df, type="prob",my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  plot_un <- assoc_plot4a(my_snp_df_un,my_corr_df, type="prob", my_locus_lab=locus, my_locus_pos=locus_pos, ylab="score", top.marker=locus)
  
  plot_name <- paste(locus, "_snps.adjusted_Model14.gassocplot.ld.index_snp.png", sep="")
  plot_name_un <- paste(locus, "_snps.unadjusted_Model14.gassocplot.ld.index_snp.png", sep="")

  assoc_plot_save(plot, plot_name)
  assoc_plot_save(plot_un, plot_name_un)

  #Create a fake correlation matrix.
  #correlation_top_snp <- matrix(NA, nrow=dim(my_snp_df)[1], ncol=dim(my_snp_df)[1]) 
  #colnames(correlation_top_snp) <- my_snp_df$marker
  #rownames(correlation_top_snp) <- my_snp_df$marker
}



