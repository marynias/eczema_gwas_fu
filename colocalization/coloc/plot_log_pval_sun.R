library(ggplot2)
library(tools)

args = commandArgs(trailingOnly=TRUE)

my_file <- args[1]
my_output <- args[2]
df <- read.table(my_file, header=T, stringsAsFactors=FALSE) 
df$log <- -log10(df$p)

fig <- ggplot(data=df, aes(df$log)) + geom_histogram(col="black") + xlab("-log10 P-value") 

ggsave(my_output, fig, dpi=300, height=5.7, width=8)