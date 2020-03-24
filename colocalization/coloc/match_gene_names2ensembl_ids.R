library(biomaRt)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}


input_file <- args[1]
output_file <- args[2]
#Read in colocalisation result as a table.
my_results <- read.delim(input_file, header=FALSE)


gene_names <- my_results$V1

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
res <- getBM(attributes = c('external_gene_name',
                            'ensembl_transcript_id', 
                            'ensembl_gene_id', 
                            'description',
                            'start_position',
                            'end_position',
                            'transcript_start',
                            'transcript_end'),
             filters = 'external_gene_name', 
             values = gene_names,
             mart = mart)

#Merge the two results tables
write.table(res, file=output_file, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)