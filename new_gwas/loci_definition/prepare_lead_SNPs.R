library(tools)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied.n", call.=FALSE)}

input_stats <- args[1]
gwas_name <- args[2]

my_gwas <- read.csv(input_stats, stringsAsFactors=F)

my_gwas <- my_gwas %>% arrange(CHR, POS)

my_gwas2 <- my_gwas %>% rename(N=N_samples) %>% rename(EFFECT_ALLELE=EA) %>% 
rename(NON_EFFECT_ALLELE=OA) %>% rename(EFFECT_ALLELE_FREQ=EAF) %>% 
rename(Z_SCORE=Z_score)
my_gwas2 <- my_gwas2 %>% dplyr::select(SNP, CHR, POS, EFFECT_ALLELE, NON_EFFECT_ALLELE, N, 
EFFECT_ALLELE_FREQ, BETA, SE, Z_SCORE, PVAL, RSID, MAF)

output_file <- paste0("leadSNPs.", gwas_name, ".txt")
write.table(my_gwas2, output_file, sep="\t", quote=F, row.names=F)

output_file2 <- paste0("leadSNPs.", gwas_name, ".index")
my_gwas3 <- my_gwas2 %>% dplyr::select(CHR, POS, RSID)
write.table(my_gwas3, output_file2, sep="\t", quote=F, col.names=F, row.names=F)

#Output 1 Mbp intervals around lead SNPs.
my_gwas4 <- my_gwas3
my_gwas4$start <- my_gwas4$POS - 500000
my_gwas4$start <- ifelse(my_gwas4$start < 1, 1, my_gwas4$start)
my_gwas4$end <- my_gwas4$POS + 500000
my_gwas4 <- my_gwas4 %>% dplyr::select(CHR, start, end, RSID)
output_file3 <- paste0("leadSNPs.", gwas_name, "_sorted_1Mbp.bed")
write.table(my_gwas4, output_file3, sep="\t", quote=F, col.names=F, row.names=F)

my_gwas5 <- my_gwas3
my_gwas5$END <- my_gwas5$POS + 1
my_gwas5 <- my_gwas5 %>% dplyr::select(CHR, POS, END, RSID)
output_file4 <- paste0("leadSNPs.", gwas_name, ".bed")
write.table(my_gwas5, output_file4, sep="\t", quote=F, col.names=F, row.names=F)
