# By default 1 million iterations are run

library(R2BGLiMS) # Load package
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least three arguments must be supplied (2 input files and Priority Pruner results).n", call.=FALSE)}

my_betas <- args[1]
my_dosage <- args[2]
pp_results <- args[3]

my_betas_df <- read.table(my_betas, header=T, stringsAsFactors=FALSE)
my_dosage_df <- read.table(my_dosage, header=T, row.names=1, stringsAsFactors=FALSE)
pp_results_df <- read.table(pp_results, header=T, stringsAsFactors=FALSE)
#Filter only SNPs with selected == 1 (ie. not pruned)
pp_results_sel <- subset(pp_results_df, selected == "1")
my_betas_df <- my_betas_df[my_betas_df$rsid %in% pp_results_sel$name, ]
my_dosage_df <- subset(my_dosage_df, rownames(my_dosage_df) %in% pp_results_sel$name)

my_betas_vector <- structure(as.numeric(my_betas_df$zscore), names = as.character(my_betas_df$rsid))
#Transpose and convert to matrix
my_dosage_matrix <- as.matrix(t(my_dosage_df))

n=103066
#n=1000
jam.results <- JAM(
  
  marginal.betas=my_betas_vector,
  
  X.ref=my_dosage_matrix,
  
  model.space.prior = list("a"=1, "b"=length(names(my_betas_vector)), "Variables"=names(my_betas_vector)),
  
  n=n)

results <- as.data.frame(jam.results@posterior.summary.table)
output_file <-  paste(file_path_sans_ext(my_betas), ".jam", sep="")
write.table(results, file=output_file, sep="\t",quote=FALSE, row.names=TRUE, col.names=TRUE)

output_pdf <-  paste(file_path_sans_ext(my_betas), ".jam.pdf", sep="")
pdf(output_pdf)
ManhattanPlot(jam.results)
dev.off()