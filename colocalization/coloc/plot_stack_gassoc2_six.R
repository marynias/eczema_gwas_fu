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
my_trait3 <- args[3]
my_trait4 <- args[4]
my_trait5 <- args[5]
my_trait6 <- args[6]
name1 = args[7]
name2 = args[8]
name3 = args[9]
name4 = args[10]
name5 = args[11]
name6 = args[12]
filename = args[13]

my_t1_df <- read.table(my_trait1, header=F, stringsAsFactors=FALSE) 
my_t2_df <- read.table(my_trait2, header=F, stringsAsFactors=FALSE) 
my_t3_df <- read.table(my_trait3, header=F, stringsAsFactors=FALSE)
my_t4_df <- read.table(my_trait4, header=F, stringsAsFactors=FALSE)
my_t5_df <- read.table(my_trait5, header=F, stringsAsFactors=FALSE)
my_t6_df <- read.table(my_trait6, header=F, stringsAsFactors=FALSE)
my_header_1 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name1)
my_t1_df  <- rename(my_t1_df, my_header_1)
my_header_2 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name2)
my_t2_df  <- rename(my_t2_df, my_header_2)
my_header_3 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name3)
my_t3_df  <- rename(my_t3_df, my_header_3)
my_header_4 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name4)
my_t4_df  <- rename(my_t4_df, my_header_4)
my_header_5 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name5)
my_t5_df  <- rename(my_t5_df, my_header_5)
my_header_6 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name6)
my_t6_df  <- rename(my_t6_df, my_header_6)

all_results2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(my_t1_df, my_t2_df, my_t3_df, my_t4_df, my_t5_df, my_t6_df))
merged_results <- all_results2[c(1, 4, 5, 6, 7, 8, 9)]
reference <-  all_results2[c(1, 2, 3)]
reference <- reference[order(reference$marker),]
markers <-  merged_results[order(merged_results$marker),]

z <- merged_results[,-1]
rownames(z) <- make.names(merged_results[,1], unique=TRUE)

#Remove missing values
#z[is.na(z)] <- 0

plot <- stack_assoc_plot(reference, z, traits=c(name1, name2, name3, name4, name5, name6), type="prob")
stack_assoc_plot_save(plot, filename, 6)

