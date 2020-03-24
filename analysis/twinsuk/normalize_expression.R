library(tools)
args = commandArgs(trailingOnly=TRUE)

my_expr <- args[1]

my_expr_df <- read.delim(my_expr, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_expr_df2 <- data.frame(my_expr_df[,-1], row.names=my_expr_df[,1])

my_all <- "twinsuk_individuals_all"
my_eczema <- "twinsuk_individuals_eczema"
my_noneczema <- "twinsuk_individuals_noneczema"

my_all_df <- read.delim(my_all, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_eczema_df <- read.delim(my_eczema, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_noneczema_df <- read.delim(my_noneczema, sep="\t", stringsAsFactors=FALSE, header=TRUE)

#Transform the expression values to standard normal
for( sl in 1:dim(my_expr_df2)[1] ) {
  mat = my_expr_df2[sl,]
  mat = t(apply(mat, 1, rank, ties.method = "average"))
  mat = qnorm(mat / (ncol(my_expr_df2)+1))
  my_expr_df2[sl,] = mat
}

#Subset to eczema, noneczema and all
eczema_intersect <- my_expr_df2[,intersect(names(my_expr_df2), my_eczema_df$id)]
noneczema_intersect <- my_expr_df2[,intersect(names(my_expr_df2), my_noneczema_df$id)]
all_intersect <- my_expr_df2[,intersect(names(my_expr_df2), my_all_df$id)]

out_eczema <- paste(file_path_sans_ext(my_expr), "_eczema_normal.rpkm", sep="")
out_noneczema <- paste(file_path_sans_ext(my_expr), "_noneczema_normal.rpkm", sep="")
out_all <- paste(file_path_sans_ext(my_expr), "_all_normal.rpkm", sep="")

#Make sure order of individuals the same throughout
write.table(t(all_intersect[ , order(names(all_intersect))]), file=out_all, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(t(noneczema_intersect[ , order(names(noneczema_intersect))]), file=out_noneczema, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(t(eczema_intersect[ , order(names(eczema_intersect))]), file=out_eczema, sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)