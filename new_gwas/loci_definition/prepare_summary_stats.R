library(tools)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied.n", call.=FALSE)}

input_stats <- args[1]
gwas_name <- args[2]

my_gwas <- read.csv(input_stats, stringsAsFactors=F)

my_gwas2 <- my_gwas %>% rename(N=N_samples) %>% rename(EFFECT_ALLELE=EA) %>% 
rename(NON_EFFECT_ALLELE=OA) %>% rename(EFFECT_ALLELE_FREQ=EAF) %>% 
rename(Z_SCORE=Z_score)
my_gwas2 <- my_gwas2 %>% dplyr::select(SNP, CHR, POS, EFFECT_ALLELE, NON_EFFECT_ALLELE, N, 
EFFECT_ALLELE_FREQ, BETA, SE, Z_SCORE, PVAL, RSID, MAF)

output_file <- paste0("results.", gwas_name, ".txt")
write.table(my_gwas2, output_file, sep="\t", quote=F, row.names=F)

#Keeping minimal number of stats for FUMA to run.
my_gwas3 <- my_gwas2 %>% dplyr::select(CHR, POS, EFFECT_ALLELE, NON_EFFECT_ALLELE, N, PVAL, RSID)

output_file2 <- paste0("results.", gwas_name, "_fuma.txt")
write.table(my_gwas3, output_file2, sep="\t", quote=F, row.names=F)