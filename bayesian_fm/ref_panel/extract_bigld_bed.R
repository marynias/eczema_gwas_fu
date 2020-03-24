library(tools)

args = commandArgs(trailingOnly=TRUE)


my_var <- args[1]
my_chrom <- args[2]
my_pos <- as.integer(args[3])
my_bigld <- args[4]


BigLDres <- read.delim(my_bigld, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_select_int <- BigLDres[BigLDres$start.bp <= my_pos & BigLDres$end.bp >= my_pos,]

file_table <- paste(my_var, "bigld", "bed", sep=".")
final_table <- cbind(my_chrom, my_select_int)
write.table(final_table[,1:4], file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

