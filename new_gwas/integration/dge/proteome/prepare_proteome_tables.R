library("stringr")
library("dplyr")

my_elias1 <- read.delim("Elias2017_TableE2.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_elias1$direction <- ifelse(my_elias1$logFC >=0, 'upregulated', 'downregulated')
my_elias1$study <- "Elias2017"
my_elias1$comparison <- "flg_knockdown_ver_normal_lse"
my_elias1$summary <- paste(my_elias1$direction, " in ", my_elias1$comparison, " (", my_elias1$study, ")", sep="")
my_elias1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_elias1 <- my_elias1[my_col_order]
write.table(my_elias1, "Elias2017_TableE2_proteome.tsv", quote=F, sep="\t", row.names=F, na="")

my_molin1 <- read.delim("Molin2015_Table3.txt", stringsAsFactors = F, header=T)
my_molin1$direction <- ifelse(my_molin1$Up.downregulation.in.CHE  == "Up", 'upregulated', 'downregulated')
my_molin1$study <- "Molin2015"
my_molin1$comparison <- "lesional_hand_eczema_ver_non-AD"
my_molin1$summary <- paste(my_molin1$direction, " in ", my_molin1$comparison, " (", my_molin1$study, ")", sep="")
my_molin1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_molin1 <- my_molin1[my_col_order]
write.table(my_molin1, "Molin2015_Table3_proteome.tsv", quote=F, sep="\t", row.names=F, na="")

my_morelli1 <- read.delim("Morelli2021_exd14276-sup-0006-tables1.txt", stringsAsFactors = F, header=T)
my_morelli1$study <- "Morelli2021"
my_morelli1$comparison <- "nonlesional_AD_ver_non-AD"
my_morelli1$summary <- paste(my_morelli1$direction, " in ", my_morelli1$comparison, " (", my_morelli1$study, ")", sep="")
my_morelli1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_morelli1 <- my_morelli1[my_col_order]
write.table(my_morelli1, "Morelli2021_exd14276-sup-0006-tables1_proteome.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2019_AD_normal_blood.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_pavel1$direction <- ifelse(my_pavel1$logFC >=0, 'upregulated', 'downregulated')
#Filter for significant DEGs:
my_pavel1 <- my_pavel1[my_pavel1$FDR < 0.05,]
my_pavel1$study <- "Pavel2020"
my_pavel1$comparison <- "blood_AD_ver_non-AD"
my_pavel1$summary <- paste(my_pavel1$direction, " in ", my_pavel1$comparison, " (", my_pavel1$study, ")", sep="")
my_pavel1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_pavel1 <- my_pavel1[my_col_order]
write.table(my_pavel1, "Pavel2019_AD_normal_blood_proteome.tsv", quote=F, sep="\t", row.names=F, na="")


my_pavel1 <- read.delim("Pavel2019_AD_lesional_normal.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_pavel1$direction <- ifelse(my_pavel1$logFC >=0, 'upregulated', 'downregulated')
#Filter for significant DEGs:
my_pavel1 <- my_pavel1[my_pavel1$FDR < 0.05,]
my_pavel1$study <- "Pavel2020"
my_pavel1$comparison <- "lesional_AD_ver_non-AD"
my_pavel1$summary <- paste(my_pavel1$direction, " in ", my_pavel1$comparison, " (", my_pavel1$study, ")", sep="")
my_pavel1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_pavel1 <- my_pavel1[my_col_order]
write.table(my_pavel1, "Pavel2019_AD_lesional_normal_proteome.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2019_AD_nonlesional_normal.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_pavel1$direction <- ifelse(my_pavel1$logFC >=0, 'upregulated', 'downregulated')
#Filter for significant DEGs:
my_pavel1 <- my_pavel1[my_pavel1$FDR < 0.05,]
my_pavel1$study <- "Pavel2020"
my_pavel1$comparison <- "nonlesional_AD_ver_non-AD"
my_pavel1$summary <- paste(my_pavel1$direction, " in ", my_pavel1$comparison, " (", my_pavel1$study, ")", sep="")
my_pavel1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_pavel1 <- my_pavel1[my_col_order]
write.table(my_pavel1, "Pavel2019_AD_nonlesional_normal_proteome.tsv", quote=F, sep="\t", row.names=F, na="")