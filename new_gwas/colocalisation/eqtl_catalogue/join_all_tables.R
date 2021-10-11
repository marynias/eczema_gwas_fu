#Join results from multiple SNPs into one table.

#Gene expression
temp = list.files(pattern="*_coloc.txt")
my_files = lapply(temp, read.delim)
final_df = do.call(rbind, my_files)
write.table(final_df, file="gxp_coloc_all_results.txt", quote=F, sep="\t", row.names=F)

#Transcript expression
temp = list.files(pattern="*_txcoloc.txt")
my_files = lapply(temp, read.delim)
final_df = do.call(rbind, my_files)
write.table(final_df, file="tx_coloc_all_results.txt", quote=F, sep="\t", row.names=F)