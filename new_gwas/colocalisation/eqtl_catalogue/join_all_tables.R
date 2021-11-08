library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least 1 arguments must be supplied.n", call.=FALSE)}

gwas_name <- args[1]

#Join results from multiple SNPs into one table.

#Gene expression
my_pattern <- paste0("*", gwas_name, "*_coloc.txt")
temp = Sys.glob(my_pattern)
my_files = lapply(temp, read.delim)
final_df = do.call(rbind, my_files)
output_file <- paste(gwas_name, "_all_coloc.txt", sep="")
write.table(final_df, file=output_file, quote=F, sep="\t", row.names=F)
