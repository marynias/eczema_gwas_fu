library(biomaRt)
library(tools)

args = commandArgs(trailingOnly=TRUE)

my_results <- args[1]
my_matrix <- read.delim(my_results, sep="\t", stringsAsFactors=FALSE, header=TRUE)

gene_ids <- my_matrix$gene

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
res <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
             filters = 'ensembl_gene_id', 
             values = gene_ids,
             mart = mart)

id <- match(my_matrix$gene, res$ensembl_gene_id)
gene_names <- res[id,c(2)]
output <- cbind(gene_names, my_matrix)

out_file_gen = paste(file_path_sans_ext(my_results), "_ensembl_names", sep="")
write.table(output, file=out_file_gen, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)