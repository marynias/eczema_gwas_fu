library(tools)

args = commandArgs(trailingOnly=TRUE)

my_ensembl <- args[1]
my_rpkm <- args[2]
my_list <- args[3]

ensembl <- read.delim(my_ensembl, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_rpkm_df <- read.delim(my_rpkm, sep=" ",stringsAsFactors=FALSE, header=TRUE)
my_order <-  read.delim(my_list, sep="\t", stringsAsFactors=FALSE, header=FALSE)
my_ids <- unlist(my_order, use.names=FALSE)
#First drop not needed individuals
my_rpkm_df2 <- my_rpkm_df[my_ids]
my_rpkm_df2$Gene <- my_rpkm_df$FID
colorder <- c("Gene", my_ids)
my_rpkm_df2 <- my_rpkm_df2[colorder]

for(i in 1:nrow(my_rpkm_df2)) 
{
    my_gene_id <- my_rpkm_df2[i,1]
    my_genes <- ensembl[ensembl$gene==my_gene_id,]
    if (nrow(my_genes) > 0)
    {
    my_rpkm_df2$Gene[i] = my_genes$gene_names
	}
	else
	{
	print(my_rpkm)
	my_rpkm_df2$Gene[i] = my_gene_id
	}
    transposed_df <- t(my_rpkm_df2[i,])
    out <- paste(my_rpkm_df2$Gene[i], "_", file_path_sans_ext(my_rpkm), ".gxp", sep="")
    write.table(transposed_df[2:length(rownames(transposed_df)),], file=out, quote=F, sep="\t", col.names = F, row.names = F)
}