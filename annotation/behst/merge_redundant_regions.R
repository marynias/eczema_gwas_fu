library(tools)

args = commandArgs(trailingOnly=TRUE)


my_bed <- args[1]
my_list <- args[2]
output <- args[3]


my_bed_df <- read.delim(my_bed, sep="\t", stringsAsFactors=FALSE, header=FALSE)
my_list_df<- read.delim(my_list, sep="\t", stringsAsFactors=FALSE, header=TRUE)

my_bed_df2 <- my_bed_df[my_bed_df$V4 %in% my_list_df$snp, ]

write.table(my_bed_df2, file=output, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)