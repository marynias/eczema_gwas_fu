#Create covariate file for MatrixEQTL

my_pheno <- "twinsuk_covariates2.txt"
my_geno <- "data.all.pca"
my_all <- "twinsuk_individuals_all"
my_eczema <- "twinsuk_individuals_eczema"
my_noneczema <- "twinsuk_individuals_noneczema"

my_pheno_df <- read.delim(my_pheno, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_pheno_df <- my_pheno_df[,c("TWIN_ID", "AGE", "ECZEMA")]
my_geno_df <- read.delim(my_geno, sep="\t", stringsAsFactors=FALSE, header=TRUE)

my_all_df <- read.delim(my_all, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_eczema_df <- read.delim(my_eczema, sep="\t", stringsAsFactors=FALSE, header=TRUE)
my_noneczema_df <- read.delim(my_noneczema, sep="\t", stringsAsFactors=FALSE, header=TRUE)

merged_df <- merge(my_pheno_df, my_geno_df, by.x="TWIN_ID", by.y="sample.id")
eczema_intersect <- intersect(merged_df$TWIN_ID, my_eczema_df$id)
noneczema_intersect <- intersect(merged_df$TWIN_ID, my_noneczema_df$id)
all_intersect <- intersect(merged_df$TWIN_ID, my_all_df$id)
eczema_covariate <- merged_df[merged_df$TWIN_ID %in% eczema_intersect, ]
eczema_covariate <- eczema_covariate[-c(3)]
noneczema_covariate <- merged_df[merged_df$TWIN_ID %in% noneczema_intersect, ]
noneczema_covariate <- noneczema_covariate[-c(3)]
all_covariate <- merged_df[merged_df$TWIN_ID %in% all_intersect, ]

#Peer Covariates
all_covariate_p <- all_covariate[order(all_covariate$TWIN_ID),]
eczema_covariate_p <- eczema_covariate[order(eczema_covariate$TWIN_ID),]
noneczema_covariate_p <- noneczema_covariate[order(noneczema_covariate$TWIN_ID),]

write.table(all_covariate_p[-c(1)], file="twinsuk_peer_covariates.all", sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(eczema_covariate_p[-c(1)], file="twinsuk_peer_covariates.eczema", sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(noneczema_covariate_p[-c(1)], file="twinsuk_peer_covariates.noneczema", sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

colnames(all_covariate)[1] <- "id"
all_covariate <- all_covariate[order(all_covariate$id),]
all_covariate <- t(all_covariate)

colnames(noneczema_covariate)[1] <- "id"
noneczema_covariate <- noneczema_covariate[order(noneczema_covariate$id),]
noneczema_covariate <- t(noneczema_covariate)

colnames(eczema_covariate)[1] <- "id"
eczema_covariate <- eczema_covariate[order(eczema_covariate$id),]
eczema_covariate <- t(eczema_covariate)

write.table(all_covariate, file="twinsuk_matrixeqtl_covariates.all", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)
write.table(eczema_covariate, file="twinsuk_matrixeqtl_covariates.eczema", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)
write.table(noneczema_covariate, file="twinsuk_matrixeqtl_covariates.noneczema", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)

write.table(all_covariate[1:3,], file="twinsuk_matrixeqtl_covariates_nopca.all", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)
write.table(eczema_covariate[1:2,], file="twinsuk_matrixeqtl_covariates_nopca.eczema", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)
write.table(noneczema_covariate[1:2,], file="twinsuk_matrixeqtl_covariates_nopca.noneczema", sep="\t",quote=FALSE, row.names=TRUE, col.names=FALSE)
