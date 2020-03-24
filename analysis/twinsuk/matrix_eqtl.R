library("MatrixEQTL")
library("tools")

args = commandArgs(trailingOnly=TRUE)

useModel = modelLINEAR
SNP_file_name = args[1]
expression_file_name = args[2]
covariates_file_name = args[3]
out_prefix = args[4]
output_file_name = "test2"

pvOutputThreshold = 1
errorCovariance = numeric()

snps = SlicedData$new();
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile( SNP_file_name )

gene = SlicedData$new();
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile( expression_file_name )

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 2000
cvrt$LoadFile( covariates_file_name )

#Filter low frequency SNPs
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

#Expresion normalization
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat)

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  output_file_name = output_file_name,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

#Results
my_results <- me$all$eqtls
out_table <- paste(out_prefix, ".matrixeqtl", sep="")
write.table(my_results, file=out_table, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)

#Pvalue plot
out_pp <- paste(out_prefix, "_pp.png", sep="")
png(filename = out_pp, width = 650, height = 650)
plot(me, col="grey")
dev.off()
# Perform the same analysis recording information for 
# a Q-Q plot
meq = Matrix_eQTL_engine(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt, 
  pvOutputThreshold = pvOutputThreshold, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  output_file_name = output_file_name,
  verbose = TRUE,
  pvalue.hist = "qqplot")

out_qq <- paste(out_prefix, "_qq.png", sep="")
png(filename = out_qq, width = 650, height = 650)
plot(meq, pch = 16, cex = 0.7)
dev.off()


