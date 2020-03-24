library(Matrix)
library("tools")
args = commandArgs(trailingOnly=TRUE)

#Results of IBD calculation
my_data <- read.delim("data.all.ibd", sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_individuals <- args[1]
my_individuals_df <- read.delim(my_individuals, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_individuals_df<- my_individuals_df[order(my_individuals_df$id),]
my_data_matrix <-  matrix(ncol = length(my_individuals_df), nrow = length(my_individuals_df))

colnames(my_data_matrix) <- my_individuals_df
rownames(my_data_matrix) <- my_individuals_df

#Fill up the empty matrix with values.
for (row in 1:nrow(my_data))
{
  id1 <- my_data[row, "ID1"]
  id2  <- my_data[row, "ID2"]
  k <- my_data[row, "kinship"]
  if (id1 %in% my_individuals_df && id2 %in% my_individuals_df)
  {
      my_data_matrix[id1,id2] = k
      my_data_matrix[id2,id1] = k
  }
}

#Need to convert to positive definite matrix.
my_data_matrix2 <- as.matrix(my_data_matrix)
my_data_matrix2[is.na(my_data_matrix2)] <- 1
my_result <- nearPD(my_data_matrix2, corr=TRUE, keepDiag=TRUE, ensureSymmetry=FALSE)
out <- as.matrix(my_result$mat)
write.table(out, file=args[2] , sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

