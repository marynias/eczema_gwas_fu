library(tools)
args = commandArgs(trailingOnly=TRUE)
library(reshape2)
my_jam <- read.delim(args[1])
my_jam <- melt(my_jam,id.vars=c("X"))
colnames(my_jam) <- c("study_id", "rsid", "PP") 
my_jam <- my_jam[!is.na(my_jam$PP),]
my_jam$PP <- as.numeric(my_jam$PP)
my_jam <- my_jam[my_jam$PP >= 0.05,]
write.table(my_jam, args[2], quote=FALSE, row.names=FALSE, sep="\t")