library(tools)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least 1 argument must be supplied.n", call.=FALSE)}

gwas_name <- args[1]

input_file <- paste0("results.", gwas_name, ".txt")

my_gwas <- read.delim(input_file, sep="\t", stringsAsFactors=F, header=T)

cols <- colnames(my_gwas)

for (my_c in cols) {
    print(my_c)
    output <- my_gwas %>% count(.data[[my_c]]) %>% arrange(desc(n))
    out_file <- paste0(my_c, "_stats_", gwas_name, ".txt")
    write.table(output, out_file, sep="\t", quote=F, row.names=F)
}
