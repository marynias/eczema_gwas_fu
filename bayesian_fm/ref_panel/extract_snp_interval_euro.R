library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least three arguments must be supplied (input file and variant name).n", call.=FALSE)}

my_tab <- args[1]
my_var <- args[2]
my_interval <- as.integer(args[3])

#my_tab <- "chr1.locus"
my_data <- read.delim(my_tab, sep="\t", stringsAsFactors=FALSE, header=TRUE)
#my_pos_line <- my_data[my_data$rsid=="rs61813875",]
my_pos_line <- my_data[my_data$RSID==my_var,]
my_pos <- my_pos_line$POS
interval_min <- my_pos - my_interval
interval_max <- my_pos + my_interval
selected <- my_data[my_data$POS>interval_min & my_data$POS<interval_max,]
my_interval <- as.integer(my_interval)
file_table <- paste(my_var, my_interval, "locus", sep=".")
write.table(selected, file=file_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)