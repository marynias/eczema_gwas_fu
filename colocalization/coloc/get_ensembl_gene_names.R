library(biomaRt)

input_file <- "gencode.v19.annotation.ensembl"
.
my_results <- read.delim(input_file, header=F)

gene_ids <- my_results$V1

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)

res <- getBM(attributes = c('ensembl_gene_id', 
                            'external_gene_name',
                            'description',
                            'start_position',
                            'end_position'),
             filters = 'ensembl_gene_id', 
             values = gene_ids,
             mart = mart)

#Merge the two results tables
colnames(my_results) <- c("gene_id", "middle")
output <- merge(my_results, res, by.x = "gene_id", by.y = "ensembl_gene_id")
write.table(output, file="gencode.v19.annotation.ensembl.names", sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)