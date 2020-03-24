########## moloc on single locus

library(moloc)
library(plyr)
library(tools)

args = commandArgs(trailingOnly=TRUE)

my_gwas <- args[1]
my_eqtl <- args[2]
my_pqtl <- args[3]

gwas <- read.table(my_gwas, header=T, stringsAsFactors=FALSE)
eqtl <- read.table(my_eqtl, header=T, stringsAsFactors=FALSE)
pqtl <- read.table(my_pqtl, header=T, stringsAsFactors=FALSE)

#Intersect SNPs in the 3 datasets to find the shared ones.
keep <- intersect(gwas$RSID, pqtl$snp)
keep2 <- intersect(keep, eqtl$rsid)
#SNPs need to be in the same order in both data frames
pqtl_keep <- subset(pqtl, snp %in% keep2)
pqtl_keep <- pqtl_keep[order(pqtl_keep$snp),] 
gwas_keep <- subset(gwas, RSID %in% keep2)
gwas_keep <- gwas_keep[order(gwas_keep$RSID),] 
eqtl_keep <- subset(eqtl, rsid %in% keep2)
eqtl_keep <- eqtl_keep[order(eqtl_keep$rsid),] 

#Rename the relevant columns.
pqtl_keep <- rename(pqtl_keep, c("snp"="SNP", "beta"="BETA", "se"="SE",
                             "p"="PVAL","effect_allele"="A1", "other_allele"="A2", "n" = "N", "effect_allele_freq"="eaf", "gene_name" = "gene"))
eqtl_keep <- rename(eqtl_keep, c("rsid"="SNP", "pval"="PVAL", "chrom"="CHR", "pos"="POS", "zscore"="z"))
gwas_keep <- rename(gwas_keep, c("SNP"="deprec","EFFECT_ALLELE"="A1", "NON_EFFECT_ALLELE"="A2", "maf"="MAF", "Z_SCORE"="z", "RSID"="SNP"))
gwas_keep$Ncases <- 0.2245 * gwas_keep$N

#Add MAF columns to pQTL and eqTL
eqtl_keep <- cbind(gwas_keep$MAF, eqtl_keep)
eqtl_keep <- rename(eqtl_keep, c("gwas_keep$MAF"="MAF"))
pqtl_keep <- cbind(gwas_keep$MAF, pqtl_keep)
pqtl_keep <- rename(pqtl_keep, c("gwas_keep$MAF"="MAF"))
#Calculate Beta and SE for eQTLgen
eqtl_keep$VARIANCE <- (1 / sqrt(2 * eqtl_keep$MAF * (1-eqtl_keep$MAF) * (14115 + eqtl_keep$z^2)))^2
eqtl_keep$BETA <- eqtl_keep$z / sqrt(2 * eqtl_keep$MAF * (1-eqtl_keep$MAF) * (14115 + eqtl_keep$z^2))
eqtl_keep$SE <- sqrt(eqtl_keep$VARIANCE)

#Add N column to eqtlgen
eqtl_keep$N <- 14115

#Run Moloc
my_input = list(gwas_keep, eqtl_keep, pqtl_keep)
my_results <- moloc_test(my_input)

my_e <- basename(my_eqtl)
output_file <- paste(file_path_sans_ext(my_e), ".pp.moloc", sep="")
output_file2 <- paste(file_path_sans_ext(my_e), ".snp.moloc", sep="")
write.table(my_results$priors_lkl_ppa, file=output_file, quote=F, sep="\t", row.names=T)
write.table(my_results$best_snp, file=output_file2, quote=F, sep="\t", row.names=T)
