library("tools")

args = commandArgs(trailingOnly=TRUE)

my_matrixeqtl_covariates <- read.delim(args[1], sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_peer_covariates <- read.delim(args[2], sep="\t", stringsAsFactors=FALSE, header=FALSE)
  
#Choose only the first 2 covariates from peer, as per plot, along with PCs and age.
my_peer_covariates_pc <- my_peer_covariates[c(1:9)]

#Choose only the first 2 covariates from peer, along with age 
my_peer_covariates_nopc <- my_peer_covariates[c(1,8:9)]

#Format to be suitable for matrixeqtl analysis. 
my_sample_names = colnames(my_matrixeqtl_covariates)
my_peer_covariates_pc <- t(my_peer_covariates_pc)
r <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
my_peer_covariates_pc2 <- cbind(r, my_peer_covariates_pc)
colnames(my_peer_covariates_pc2) <- my_sample_names


my_peer_covariates_nopc <- t(my_peer_covariates_nopc)
r2 <- c("V1", "V2", "V3")
my_peer_covariates_nopc2 <- cbind(r2, my_peer_covariates_nopc)
colnames(my_peer_covariates_nopc2) <- my_sample_names

out_table1 <- paste(file_path_sans_ext(args[2]), "_pc_covariates", sep="")
out_table2 <- paste(file_path_sans_ext(args[2]), "_nopc_covariates", sep="")
write.table(my_peer_covariates_pc2, file=out_table1, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(my_peer_covariates_nopc2, file=out_table2, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)