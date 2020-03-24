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
name1=args[3]
name2=args[4]
filename=args[5]

my_t1_df<- read.table(my_trait1, header=F, stringsAsFactors=FALSE) 
my_t2_df<- read.table(my_trait2, header=F, stringsAsFactors=FALSE) 
my_header_1 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name1)
my_t1_df  <- rename(my_t1_df, my_header_1)
my_header_2 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"=name2)
my_t2_df  <- rename(my_t2_df, my_header_2)

all_results <- merge(my_t1_df, my_t2_df, all=TRUE)
merged_results <- all_results[c(1, 4, 5)]

reference <-  all_results[c(1, 2, 3)]
reference <- reference[order(reference$marker),]
markers <-  merged_results[order(merged_results$marker),]

z <- merged_results[,-1]
rownames(z) <- make.names(merged_results[,1], unique=TRUE)

#Remove missing values
#z[is.na(z)] <- 0

plot <- stack_assoc_plot(reference, z, traits=c(name1, name2), type="prob")
stack_assoc_plot_save(plot, filename, 2)

