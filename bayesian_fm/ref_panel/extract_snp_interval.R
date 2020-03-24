library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least three arguments must be supplied (input file and variant name).n", call.=FALSE)}

my_tab <- args[1]
my_var <- args[2]
my_interval <- as.integer(args[3])

#my_tab <- "chr1.locus"
my_data <- read.delim(my_tab, sep=" ", stringsAsFactors=FALSE)
#my_pos_line <- my_data[my_data$rsid=="rs61813875",]
my_pos_line <- my_data[my_data$rsid==my_var,]
my_pos <- my_pos_line$pos
interval_min <- my_pos - my_interval
interval_max <- my_pos + my_interval
selected <- my_data[my_data$pos>interval_min & my_data$pos<interval_max,]
file_table <- paste(file_path_sans_ext(args[1]), my_var, my_interval, "locus2", sep=".")
write.table(selected, file=file_table, sep=" ",quote=FALSE, row.names=FALSE, col.names=TRUE)