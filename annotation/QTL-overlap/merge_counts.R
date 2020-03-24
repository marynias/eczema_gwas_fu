library(reshape)
library(plyr)
library(tools)

args = commandArgs(trailingOnly=TRUE)
output <- args[1]

filenames <- list.files(pattern = "\\.counts$", full.names=TRUE)
import.list <- llply(filenames, read.delim, check.names=FALSE)


data <- merge_recurse(import.list)
summed <- rowSums(data[-1])
data$sum <- summed

write.table(data, output, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)