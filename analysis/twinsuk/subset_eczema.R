library(tools)
args = commandArgs(trailingOnly=TRUE)

my_rpkm <- args[1]
my_eczema <-args[2]

#my_rpkm <-"rs12188917_3Mbp_snps_headers.rpkm"
#my_eczema <- "chr2_Eurobats_Public_alt_pheno.sample"

my_rpkm_df <- read.delim(my_rpkm, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)
my_eczema_df <- read.delim(my_eczema, sep="\t", stringsAsFactors=FALSE, header=TRUE)
eczema = my_eczema_df[!is.na(my_eczema_df$Phenotype_1) &  my_eczema_df$Phenotype_1==1,]
noneczema = my_eczema_df[!is.na(my_eczema_df$Phenotype_1) &  my_eczema_df$Phenotype_1==0,]

eczema_expression <- my_rpkm_df[, intersect(names(my_rpkm_df), eczema$ID_1)] 
noneczema_expression <- my_rpkm_df[, intersect(names(my_rpkm_df), noneczema$ID_1)] 

out_eczema <- paste(file_path_sans_ext(my_rpkm), "_eczema.rpkm", sep="")
out_noneczema <- paste(file_path_sans_ext(my_rpkm), "_noneczema.rpkm", sep="")
out_all <- paste(file_path_sans_ext(my_rpkm), "_all.rpkm", sep="")
write.table(eczema_expression, file=out_eczema, sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)
write.table(noneczema_expression, file=out_noneczema, sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)
write.table(my_rpkm_df, file=out_all, sep="\t",quote=FALSE, row.names=TRUE, col.names=NA)