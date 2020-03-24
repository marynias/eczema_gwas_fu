library(coloc)

#generate input files for coloc analysis.

#read in results from eQTL analysis 

eQTL_data <- read.table("CD207_eQTL_analysis_test.assoc_sorted.txt", header=T)

#define region around lead SNP (100kb either side of SNP)
SNP_pos <- eQTL_data$ps[which(eQTL_data$rs == "rs112111458")]
upp_pos <- SNP_pos + 100000
low_pos <- SNP_pos - 100000


#read in GWAS data and subset for SNPs from eQTL analysis
#plus 23andme estimates
ECZ_GWAS <- read.table("/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015/results.euro.tsv", header=T)

#calculate z-scores and p-values:
ECZ_GWAS$Z_SCORE <- ECZ_GWAS$BETA/ECZ_GWAS$SE
ECZ_GWAS$P_VALUE <- 2*pnorm(-abs(ECZ_GWAS$Z_SCORE))

#combine GWAs and eQTL data
colnames(eQTL_data)[3] <- "POS"
GWAS_data <- merge(ECZ_GWAS, eQTL_data, by="POS")

locus_snps <- GWAS_data[which(GWAS_data$POS > low_pos & GWAS_data$POS < upp_pos), ]


#Perform coloc


#N = number of people (102762 for GWAS dataset including 23andme, 672 for eQTL dataset)
#N minus 23andme = 102762 - 62231 = 40531
#ratio of cases to controls = 0.183
#using p-values:

my.res <- coloc.abf(dataset1=list(pvalues=locus_snps$P_VALUE, N=102762, snp=locus_snps$rs, type="cc", s=0.183),
            dataset2=list(pvalues=locus_snps$p_score, N=672, snp=locus_snps$rs, type="quant"),
            MAF=locus_snps$af)

#view summary
my.res$summary

#order by decreasing p-value and save results to file
coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]

write.table(coloc_results, file="results_CD207_plus_23andme.txt-", quote=F, sep=" ", row.names=F)