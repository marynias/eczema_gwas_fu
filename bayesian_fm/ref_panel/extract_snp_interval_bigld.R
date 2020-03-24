library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
  stop("At least five arguments must be supplied (input file and variant name).n", call.=FALSE)}


my_var <- args[1]
my_chrom <- args[2]
my_pos <- as.integer(args[3])
my_bigld <- args[4]
my_tab <- args[5]

BigLDres <- read.delim(my_bigld, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_select_int <- BigLDres[BigLDres$start.bp <= my_pos & BigLDres$end.bp >= my_pos,]


my_data <- read.delim(my_tab, sep="\t", stringsAsFactors=FALSE, header=TRUE)
selected <- my_data[my_data$POS>=my_select_int$start.bp[1] & my_data$POS<=my_select_int$end.bp[1] ,]
selected <- selected[selected$CHR==my_chrom,]
file_table <- paste(my_var, "locus", sep=".")
write.table(selected, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)