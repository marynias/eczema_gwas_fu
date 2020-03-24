library(dplyr)
library(tidyr)
#Read in sample SNP table.
sample_snp <- read.delim("rs112111458_rs112111458.snps", header=T, sep="\t")
sample_snp <- tbl_df(sample_snp)
#Read in sample gene table.
sample_gene <- read.delim("rs112111458_CD207.genes", header=T, sep="\t")
sample_gene <- tbl_df(sample_gene)


#Read in locus names
locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t")
colnames(locus_names)[2] <- "locus"

#Read in table descriptions
table_desc <- read.delim("Resource_index.tsv", header=T, sep="\t")



#Read in top ranked genes for each locus - Model13.
model13_gene <- read.delim("Model13_gene_ranked.txt", sep="\t", header=T)
model13_gene_top3 <- model13_gene[model13_gene$rank < 4,]
#Read in top ranked SNP for each locus - Model09.
model09_snp <- read.delim("Model09_snp_ranked.txt", sep="\t", header=T)
model09_snp_top3 <- model09_snp[model09_snp$rank < 4,]

#Read in top ranked SNP for each locus - Model10.
model10_snp <- read.delim("Model10_snp_ranked.txt", sep="\t", header=T)
model10_snp_top3 <- model10_snp[model10_snp$rank < 4,]

model14_snp <- read.delim("Model14_snp_ranked.txt", sep="\t", header=T)
model14_snp_top3 <- model14_snp[model14_snp$rank < 4,]


#Read in the relevant results
model09_snp_top3$file <- paste("Model09_top_10_snp/", model09_snp_top3$index_SNP_rsid, '_', model09_snp_top3$current_SNP_rsid, ".snps", sep="")

model10_snp_top3$file <- paste("Model10_top_10_snp/", model10_snp_top3$index_SNP_rsid, '_', model10_snp_top3$current_SNP_rsid, ".snps", sep="")

model14_snp_top3$file <- paste("Model14_top_10_snp/", model14_snp_top3$index_SNP_rsid, '_', model14_snp_top3$current_SNP_rsid, ".snps", sep="")

model13_gene_top3$file <- paste("Model13_top_10_gene/", model13_gene_top3$index_SNP_rsid, '_', model13_gene_top3$gene_name, ".genes", sep="")

model_all_snps <- rbind(model09_snp_top3, model10_snp_top3, model14_snp_top3)


##SNPs
for (my_file in unique(model_all_snps$file))
{
  sample_snp <- read.delim(my_file, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  sample_snp <- tbl_df(sample_snp)
  #Mutate study_names
  sample_snp <- sample_snp %>% mutate(study_type = recode(study_type, 'oe' = 'enhancer', 'promoter enhancer'='promoter-enhancer','y'='enhancer','enhancers'='enhancer','bait'='promoter', 'pir'='enhancer','x'='promoter','regulatory variant predictions'='regulatory variant prediction', 'enhancer-promoter'='promoter-enhancer', 'coloc - SNP'='coloc', 'FAIRE'='FAIRE-seq', 'DEG'='DGE', 'CTCF sites'='CTCF', 'coloc - gene'='coloc', 'network analysis'='PrixFixe', 'ChIP-seq peaks'='ChIP-Seq', 'epigenome'='active chromatin state', 'Converved genomic regulatory blocks'='Conserved genomic regulatory blocks', 'ChIP-Seq'='ChIP-seq'))
  
  #Drop unnecessary columns
  sample_snp <- sample_snp[c("table_id","study_id", "index_SNP_rsid", "current_SNP_rsid", "gene_name", "beta", "FDR", "p.value", "posterior_prob", "score", "effect_allele", "tissue", "sample_size", "study_type", "evidence_weight")]
  
  sample_snp <-  merge(sample_snp, table_desc[c("table_ID", "table_description")], by.x="table_id", by.y="table_ID", all.x=TRUE)
  sample_snp <- sample_snp[order(sample_snp$evidence_weight),]
  sample_snp <- sample_snp %>% unite("FDR/p.value/posterior_probability/score", c("FDR", "p.value", "posterior_prob", "score"), sep="/", remove=FALSE)
  sample_snp <- sample_snp[-c(8:11)]
  sample_snp <- merge(sample_snp, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=TRUE)
  sample_snp <- sample_snp %>% select("table_id", "study_id", "table_description", "index_SNP_rsid", "locus", "current_SNP_rsid", "gene_name", "FDR/p.value/posterior_probability/score", "beta", "effect_allele", "tissue", "sample_size", "study_type", "evidence_weight")
  my_output <- paste(my_file, ".summary", sep="")
  write.table(sample_snp, my_output, sep="\t", col.names=T, row.names=F, quote=F)
}


for (my_file in unique(model13_gene_top3$file))
{
  sample_gene <- read.delim(my_file, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  sample_gene <- tbl_df(sample_gene)
  sample_gene <- sample_gene %>% mutate(study_type = recode(study_type, 'oe' = 'enhancer', 'promoter enhancer'='promoter-enhancer','y'='enhancer','enhancers'='enhancer','bait'='promoter', 'pir'='enhancer','x'='promoter','regulatory variant predictions'='regulatory variant prediction', 'enhancer-promoter'='promoter-enhancer', 'coloc - SNP'='coloc', 'FAIRE'='FAIRE-seq', 'DEG'='DGE', 'CTCF sites'='CTCF', 'coloc - gene'='coloc', 'network analysis'='PrixFixe', 'ChIP-seq peaks'='ChIP-Seq', 'epigenome'='active chromatin state', 'Converved genomic regulatory blocks'='Conserved genomic regulatory blocks', 'ChIP-Seq'='ChIP-seq'))
  sample_gene <- sample_gene[c("table_id","study_id", "index_SNP_rsid", "current_SNP_rsid", "gene_name", "beta", "FDR", "p.value", "posterior_prob", "score", "effect_allele", "tissue", "sample_size", "study_type","evidence_weight")]
  #Pick best p-value/FDR for a given table.
  top_fdr <- sample_gene %>% group_by(table_id, tissue) %>% top_n(-1, FDR) %>% ungroup
  top_pval <- sample_gene %>% group_by(table_id, tissue) %>% top_n(-1, p.value) %>% ungroup
  #For overlaps based on position.
  other <- sample_gene %>% filter(is.na(FDR)) %>% filter(is.na(p.value)) %>% filter(is.na(posterior_prob)) %>% filter(is.na(score))
  #Get one representative hit.
  other <- other %>% group_by(table_id, tissue) %>% top_n(1, current_SNP_rsid) %>% ungroup
  #Concatenate the two, filter for unique
  concat <- unique(rbind(top_fdr, top_pval))
  #The complement
  multiple_hits <- rbind(other, concat)
  complement <- subset(sample_gene, !(table_id %in% multiple_hits$table_id))
  altogether <- rbind(multiple_hits, complement)
  #Add counts of hits.
  counts <- sample_gene %>% group_by(table_id, tissue) %>% summarise(count = n())
  merged_counts <- merge(altogether, counts[c("table_id", "count")], by="table_id", all.x=TRUE)
  
  #Add table descriptions
  merged_counts_descriptions <- merge(merged_counts, table_desc[c("table_ID", "table_description")], by.x="table_id", by.y="table_ID", all.x=TRUE)
  merged_counts_descriptions <- merged_counts_descriptions[order(merged_counts_descriptions$evidence_weight),]
  merged_columns <- merged_counts_descriptions %>% unite("FDR/p.value/posterior_probability/score", c("FDR", "p.value", "posterior_prob", "score"), sep="/", remove=FALSE)
  merged_columns <- merged_columns[-c(9:12)]
  merged_columns <- merge(merged_columns, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=TRUE)
  merged_columns <- merged_columns %>% select("table_id", "study_id", "table_description", "index_SNP_rsid", "locus", "current_SNP_rsid", "gene_name", "FDR/p.value/posterior_probability/score", "beta", "effect_allele", "tissue", "sample_size", "study_type", "evidence_weight", "count")
  my_output <- paste(my_file, ".summary", sep="")
  write.table(merged_columns, my_output, sep="\t", col.names=T, row.names=F, quote=F)
}