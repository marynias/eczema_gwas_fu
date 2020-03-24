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
filename <- args[2]

my_t1_df<- read.table(my_trait1, header=F, stringsAsFactors=FALSE) 

my_header_1 <- c("V1"="marker", "V2"="chr", "V3"="pos", "V4"="prob")
my_t1_df <- rename(my_t1_df, my_header_1)

my_reference <- my_t1_df[c(1,2,3)]
my_z <- my_t1_df[,-1]
rownames(my_z) <- make.names(my_t1_df[,1], unique=TRUE)
my_z <- my_z[c(3)]

plot <- stack_assoc_plot(my_reference, my_z, traits=c("eQTL"), type="prob")
stack_assoc_plot_save(plot, filename, 1)