library(gassocplot)
library(gtable)
library(ggplot2)
library(grid) 
library(gridExtra)
library(plyr)
library(tools)

args = commandArgs(trailingOnly=TRUE)

### stack_assoc_plot
my_trait1 <- args[1]
my_trait2 <- args[2]
my_gwas <- args[3]
name1=args[4]
name2=args[5]
filename=args[6]

my_t1_df<- read.table(my_trait1, header=F, stringsAsFactors=FALSE) 
my_t2_df<- read.table(my_trait2, header=F, stringsAsFactors=FALSE) 
my_gwas_df <- read.table(my_gwas, header=T, stringsAsFactors=FALSE) 

#Keep only the top hit with biggest probability.
my_t1_df_max <- my_t1_df[which.max(my_t1_df$V4),]
my_header_1 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name1)
my_t1_df_max  <- rename(my_t1_df_max , my_header_1)
my_t2_df_max <- my_t2_df[which.max(my_t2_df$V4),]
my_header_2 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name2)
my_t2_df_max  <- rename(my_t2_df_max , my_header_2)

merged_results <- merge(my_t1_df_max, my_t2_df_max, all=TRUE)
merged_results <- merged_results[c(1, 4, 5)]

my_gwas_df <- read.table(my_gwas, header=T, stringsAsFactors=FALSE) 
#Drop not needed columns
my_gwas_df <- my_gwas_df[c(12, 2, 3)]
my_header_gwas <- c("RSID"="marker", "CHR"="chr", "POS"="pos")
my_gwas_df  <- rename(my_gwas_df , my_header_gwas)
markers <- my_gwas_df[order(my_gwas_df$marker),]

final_merge <- merge(my_gwas_df, merged_results, by="marker",all.x=TRUE)
final_merge[is.na(final_merge)] <- 0
final_merge <-  final_merge[c(1, 4, 5)]
final_merge <- final_merge[order(final_merge$marker),]

z <- final_merge[,-1]
rownames(z) <- final_merge[,1]

plot <- stack_assoc_plot(markers, z, traits=c(name1, name2), type="prob")
stack_assoc_plot_save(plot, filename, 2)