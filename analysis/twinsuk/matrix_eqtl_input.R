library("tools")

args = commandArgs(trailingOnly=TRUE)

my_snp <-args[1]
my_sample <- args[2]
my_expr <- args[3]
my_all <- args[4]

my_snp_df <- read.delim(my_snp, sep=" ", stringsAsFactors=FALSE, header=FALSE)
my_sample_df <- read.delim(my_sample, sep=" ", stringsAsFactors=FALSE, header=TRUE)
my_expr_df <- read.delim(my_expr, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_expr_df2 <- data.frame(my_expr_df[,-1], row.names=my_expr_df[,1])
my_all_df <- read.delim(my_all, sep="\t", stringsAsFactors=FALSE, header=TRUE)

#Subset to samples with genotypes 
my_expr_df2 <- my_expr_df2[,intersect(colnames(my_expr_df2), my_all_df$id)]

#Ignore first row
my_sample_df <- my_sample_df[2:dim(my_sample_df)[1],]
sample_names <- my_sample_df$ID_1

gene_names <- rownames(my_expr_df2)
rownames(my_expr_df2) <- c()
#output expression
my_expr_df3 <- cbind(id = gene_names, my_expr_df2) 
out_expr <- paste(file_path_sans_ext(my_snp), ".expression", sep="")
write.table(my_expr_df3[ , order(names(my_expr_df3))], file=out_expr, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)

#Output genotypes.
my_row_names <- my_snp_df[,c(2)]
my_snp_df <- my_snp_df[,c(6:dim(my_snp_df)[2])]
colnames(my_snp_df) <- sample_names

#Filter genotypes to only output those where we have expression.
to_keep <- my_all_df$id

my_snp_df<-my_snp_df[,intersect(to_keep, colnames(my_snp_df))]
my_snp_df2 <- cbind(id=my_row_names, my_snp_df)

out_geno <- paste(file_path_sans_ext(my_snp), ".genotype", sep="")
write.table(my_snp_df2[ , order(names(my_snp_df2))], file=out_geno, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)