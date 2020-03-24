library(biomaRt)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}


input_file <- args[1]
#Read in colocalisation result as a table.
my_results <- read.delim(input_file)
transcript_ids <- my_results$Ensemble_ID

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
res <- getBM(attributes = c('ensembl_transcript_id', 
                            'ensembl_gene_id', 
                            'external_gene_name',
                            'description',
                            'start_position',
                            'end_position',
                            'transcript_start',
                            'transcript_end'),
             filters = 'ensembl_transcript_id', 
             values = transcript_ids,
             mart = mart)

#Merge the two results tables
output <- merge(my_results, res, by.x = "Ensemble_ID", by.y = "ensembl_transcript_id")
colnames(output)[1] <- "ensembl_transcript_id"
#Write full results, with all affected isoforms
out_file_trans <- paste(file_path_sans_ext(input_file), ".ensembl.transcript", sep="")
write.table(output, file=out_file_trans, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)
#Write filtered results, with only affected genes (for a given SNP).
output_genes <- output[!duplicated(output[c(2,7)]), ]
out_file_gene <- paste(file_path_sans_ext(input_file), ".ensembl.gene", sep="")
write.table(output_genes, file=out_file_gene, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)