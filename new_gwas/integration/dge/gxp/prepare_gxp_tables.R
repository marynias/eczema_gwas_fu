library("stringr")
library("dplyr")

my_cole1 <- read.delim("Table_E1_Cole_(2014)_Filaggrin-stratified_transcriptomic_analysis_supp1.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_cole1$direction <- ifelse(my_cole1$logFC >=0, 'upregulated', 'downregulated')
#Filter for significant DEGs:
my_cole1 <- my_cole1[my_cole1$FDR < 0.05,]
my_cole1$summary <- paste(my_cole1$direction, " in ", my_cole1$comparison, " (", my_cole1$study, ")", sep="")
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_cole1 <- my_cole1[my_col_order]
write.table(my_cole1, "Cole2014_TableE1_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_ewald1 <- read.delim("Ewald2015_12920_2015_133_MOESM7_ESM.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_ewald1$direction <- ifelse(my_ewald1$logFC >=0, 'upregulated', 'downregulated')
#Filter for significant DEGs:
my_ewald1 <- my_ewald1[my_ewald1$FDR < 0.05,]
my_ewald1$summary <- paste(my_ewald1$direction, " in ", my_ewald1$comparison, " (", my_ewald1$study, ")", sep="")
my_ewald1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_ewald1 <- my_ewald1[my_col_order]
write.table(my_ewald1, "Ewald2015_TableM7_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_winge1 <- read.delim("Winge2011_Supp_Tab_S1_flg.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_winge1$direction <- ifelse(my_winge1$logFC >=0, 'upregulated', 'downregulated')
my_winge1$study <- "Winge2011"
my_winge1$comparison <- "nonlesional_AD_flg_homozygote_ver_non-AD"
my_winge1$summary <- paste(my_winge1$direction, " in ", my_winge1$comparison, " (", my_winge1$study, ")", sep="")
my_winge1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_winge1 <- my_winge1[my_col_order]
write.table(my_winge1, "Winge2011_Supp_Tab_S1_flg_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_winge1 <- read.delim("Winge2011_Supp_Tab_S1_flg_het.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_winge1$direction <- ifelse(my_winge1$logFC >=0, 'upregulated', 'downregulated')
my_winge1$study <- "Winge2011"
my_winge1$comparison <- "nonlesional_AD_flg_heterozygote_ver_non-AD"
my_winge1$summary <- paste(my_winge1$direction, " in ", my_winge1$comparison, " (", my_winge1$study, ")", sep="")
my_winge1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_winge1 <- my_winge1[my_col_order]
write.table(my_winge1, "Winge2011_Supp_Tab_S1_flg_het_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_winge1 <- read.delim("Winge2011_Supp_Tab_S1_flg_no.txt", stringsAsFactors = F, header=T)
#Add direction of expression change
my_winge1$direction <- ifelse(my_winge1$logFC >=0, 'upregulated', 'downregulated')
my_winge1$study <- "Winge2011"
my_winge1$comparison <- "nonlesional_AD_flg_wt_ver_non-AD"
my_winge1$summary <- paste(my_winge1$direction, " in ", my_winge1$comparison, " (", my_winge1$study, ")", sep="")
my_winge1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_winge1 <- my_winge1[my_col_order]
write.table(my_winge1, "Winge2011_Supp_Tab_S1_flg_no_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_tsoi1 <- read.delim("Tsoi2020_acute_chronic_AD_lesions.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_tsoi1 <- my_tsoi1[my_tsoi1$FDR < 0.05,]
#Add direction of expression change
my_tsoi1$direction <- ifelse(my_tsoi1$logFC >=0, 'upregulated', 'downregulated')
my_tsoi1$study <- "Tsoi2020"
my_tsoi1$comparison <- "lesional_acute_ver_chronic_AD"
my_tsoi1$summary <- paste(my_tsoi1$direction, " in ", my_tsoi1$comparison, " (", my_tsoi1$study, ")", sep="")
my_tsoi1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_tsoi1 <- my_tsoi1[my_col_order]
write.table(my_tsoi1, "Tsoi2020_acute_chronic_AD_lesions_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_rohajn1 <-  read.delim("Rojahn2020_TableE8_AD_control_blisters_upregulated.txt", stringsAsFactors = F, header=T)
my_rohajn2 <-  read.delim("Rojahn2020_TableE9_control_blisters_downregulated.txt", stringsAsFactors = F, header=T)
my_rojahn_combined <- rbind(my_rohajn1, my_rohajn2)
#Filter for significant DEGs:
my_rojahn_combined <- my_rojahn_combined[my_rojahn_combined$FDR < 0.05,]
#Add direction of expression change
my_rojahn_combined$direction <- ifelse(my_rojahn_combined$logFC >=0, 'upregulated', 'downregulated')
my_rojahn_combined$study <- "Rojahn2020"
my_rojahn_combined$comparison <- "lesional_AD_ver_non-AD"
my_rojahn_combined$summary <- paste(my_rojahn_combined$direction, " in ", my_rojahn_combined$comparison, " (", my_rojahn_combined$study, ")", sep="")
my_rojahn_combined$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_rojahn_combined <- my_rojahn_combined[my_col_order]
write.table(my_rojahn_combined, "Rojahn2020_lesional_AD_ver_non-AD_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2020_Table_S3_lesional_vs_normal.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_pavel1 <- my_pavel1[my_pavel1$FDR < 0.05,]
#Add direction of expression change
my_pavel1$direction <- ifelse(my_pavel1$logFC >=0, 'upregulated', 'downregulated')
my_pavel1$study <- "Pavel2021"
my_pavel1$comparison <- "paedriatic_lesional_AD_ver_non-AD"
my_pavel1$summary <- paste(my_pavel1$direction, " in ", my_pavel1$comparison, " (", my_pavel1$study, ")", sep="")
my_pavel1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_pavel1 <- my_pavel1[my_col_order]
write.table(my_pavel1, "Pavel2020_Table_S3_lesional_vs_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2020_Table_S3_nonlesional_vs_normal.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_pavel1 <- my_pavel1[my_pavel1$FDR < 0.05,]
#Add direction of expression change
my_pavel1$direction <- ifelse(my_pavel1$logFC >=0, 'upregulated', 'downregulated')
my_pavel1$study <- "Pavel2021"
my_pavel1$comparison <- "paedriatic_nonlesional_AD_ver_non-AD"
my_pavel1$summary <- paste(my_pavel1$direction, " in ", my_pavel1$comparison, " (", my_pavel1$study, ")", sep="")
my_pavel1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_pavel1 <- my_pavel1[my_col_order]
write.table(my_pavel1, "Pavel2020_Table_S3_nonlesional_vs_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_he1 <- read.delim("He2021_lesional_AD_ver_normal.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_he1 <- my_he1[my_he1$FDR < 0.05,]
#Add direction of expression change
my_he1$direction <- ifelse(my_he1$logFC >=0, 'upregulated', 'downregulated')
my_he1$study <- "He2021"
my_he1$comparison <- "lesional_AD_ver_non-AD"
my_he1$summary <- paste(my_he1$direction, " in ", my_he1$comparison, " (", my_he1$study, ")", sep="")
my_he1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_he1 <- my_he1[my_col_order]
write.table(my_he1, "He2021_lesional_AD_ver_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_he1 <- read.delim("He2021_nonlesional_AD_ver_normal.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_he1 <- my_he1[my_he1$FDR < 0.05,]
#Add direction of expression change
my_he1$direction <- ifelse(my_he1$logFC >=0, 'upregulated', 'downregulated')
my_he1$study <- "He2021"
my_he1$comparison <- "nonlesional_AD_ver_non-AD"
my_he1$summary <- paste(my_he1$direction, " in ", my_he1$comparison, " (", my_he1$study, ")", sep="")
my_he1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_he1 <- my_he1[my_col_order]
write.table(my_he1, "He2021_nonlesional_AD_ver_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_dyjack1 <- read.delim("Dyjack2018_TableE4_ADlesional_ver_healthy_control.txt", stringsAsFactors = F, header=T)
#Filter for significant DEGs:
my_dyjack1 <- my_dyjack1[my_dyjack1$FDR < 0.05,]
#Add direction of expression change
my_dyjack1$direction <- ifelse(my_dyjack1$logFC >=0, 'upregulated', 'downregulated')
my_dyjack1$study <- "Dyjack2018"
my_dyjack1$comparison <- "lesional_AD_ver_non-AD"
my_dyjack1$summary <- paste(my_dyjack1$direction, " in ", my_dyjack1$comparison, " (", my_dyjack1$study, ")", sep="")
my_dyjack1$Ensembl_gene_ID <- NA
my_col_order <- c("Ensembl_gene_ID",	"HGNC_symbol",	"summary")
my_dyjack1 <- my_dyjack1[my_col_order]
write.table(my_dyjack1, "Dyjack2018_TableE4_ADlesional_ver_healthy_control_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2019_AD_lesional_normal_rnaseq.txt", stringsAsFactors = F, header=T)
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
write.table(my_pavel1, "Pavel2019_AD_lesional_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")

my_pavel1 <- read.delim("Pavel2019_AD_nonlesional_normal_rnaseq.txt", stringsAsFactors = F, header=T)
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
write.table(my_pavel1, "Pavel2019_AD_nonlesional_normal_gxp.tsv", quote=F, sep="\t", row.names=F, na="")
