library(tools)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least 1 argument must be supplied.n", call.=FALSE)}

gwas_name <- args[1]

search_pattern <- paste0(gwas_name, "_open_targets_part*", "postgap.txt")
all_result_files <- Sys.glob(search_pattern)

loadFile <- function(x) {
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  print(length(colnames(df)))
  return(df)
}

all_normal <- lapply(all_result_files, loadFile)
all_normal_together <- do.call(rbind,all_normal)
all_normal_together <- as.data.frame(all_normal_together)

output <- paste0(gwas_name, "_postgap.txt")
write.table(all_normal_together, output, sep="\t", quote=F, row.names=F, col.names=T)