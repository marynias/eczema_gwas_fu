library("biomaRt")
expr_table <- read.table("GSM815426_ANL.txt", header=T, stringsAsFactors=FALSE, comment.char="#")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart=mart, attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "gene_biotype", "hgnc_symbol"), filter="affy_hg_u133_plus_2", values=expr_table$ID_REF, uniqueRows=TRUE)
#annotLookup <- read.table("affy_ids", header=T, stringsAsFactors=FALSE, fill=TRUE)
annotLookup$hgnc_symbol <- ifelse(annotLookup$hgnc_symbol == "", annotLookup$ensembl_gene_id, annotLookup$hgnc_symbol)

normal_files <- list.files(pattern="*_Normal.txt$", recursive=F)
lesion_files <-list.files(pattern="*_AL.txt$", recursive=F)
nonlesion_files <- list.files(pattern="*_ANL.txt$", recursive=F)

loadFile <- function(x) {
df <- read.delim(x, header=T, stringsAsFactors=FALSE, row.names=1, comment.char="#")
sample_name <-sub("(\\w+)_\\w+\\.txt", "\\1", basename(x))
colnames(df) <- c(sample_name)
df
}

all_normal <- lapply(normal_files, loadFile)
all_normal_master <- do.call(cbind, all_normal)

all_lesion <- lapply(lesion_files, loadFile)
all_lesion_master <- do.call(cbind, all_lesion)

all_nonlesion <- lapply(nonlesion_files, loadFile)
all_nonlesion_master <- do.call(cbind, all_nonlesion)

#Substitute Affymetrix IDs with gene names and sum gene expression values for probes matching the same gene.
process_df <- function(x)
{
gene_name <- annotLookup$hgnc_symbol[match(rownames(x), annotLookup$affy_hg_u133_plus_2)]
x$gene_name <- as.factor(gene_name)
my_temp <- lapply(x[-c(-1)], as.numeric)
x[-c(-1)] <- do.call(cbind, my_temp)
x <- aggregate(. ~ gene_name, x, sum)
}
all_normal_master <- process_df(all_normal_master)
all_lesion_master <- process_df(all_lesion_master)
all_nonlesion_master <- process_df(all_nonlesion_master)

write.table(all_normal_master, file="all_normal.expr", quote=F, sep="\t", col.names=T, row.names=F)
write.table(all_lesion_master, file="all_lesion.expr", quote=F, sep="\t", col.names=T, row.names=F)
write.table(all_nonlesion_master, file="all_nonlesion.expr", quote=F, sep="\t", col.names=T, row.names=F)