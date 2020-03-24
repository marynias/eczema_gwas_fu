library(gassocplot)
library(gtable)
library(ggplot2)
library(grid) 
library(gridExtra)
library(plyr)
library(tools)

args = commandArgs(trailingOnly=TRUE)


###Single plot

my_trait1 <- args[1]
my_gwas <- args[2]
filename <- args[3]

my_t1_df<- read.table(my_trait1, header=F, stringsAsFactors=FALSE) 
my_gwas_df <- read.table(my_gwas, header=T, stringsAsFactors=FALSE) 
#Drop not needed columns
my_gwas_df <- my_gwas_df[c(12, 2, 3)]
my_header_gwas <- c("RSID"="marker", "CHR"="chr", "POS"="pos")
my_gwas_df  <- rename(my_gwas_df , my_header_gwas)#Keep only the top hit with biggest probability.

#Keep only top pvalue
my_t1_df_max <- my_t1_df[which.max(my_t1_df$V4),]
my_header_1 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"="prob")
my_t1_df_max  <- rename(my_t1_df_max , my_header_1)

merged_gwas <- merge(my_gwas_df, my_t1_df_max, by="marker", all.x=TRUE)
colnames(merged_gwas)[2] = "chr"
colnames(merged_gwas)[3] = "pos"
merged_gwas <- merged_gwas[c(1, 2, 3, 6)]
#Convert missing values to 0 - only way to produce a plot
merged_gwas[is.na(merged_gwas)] <- 0
markers <- merged_gwas[order(merged_gwas$marker),]
plot <- assoc_plot(markers, type="prob")
assoc_plot_save(plot, filename)