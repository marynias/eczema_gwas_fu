library("dplyr")
library("stringr")
my_gene_ref <- read.table("../paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.hugo", header=T, stringsAsFactors=TRUE, fill=TRUE)
my_snp_ref <- read.table("../interval_r2_0.2_1k_nodups_sync.map", header=F, stringsAsFactors=TRUE, fill=TRUE)
colnames(my_snp_ref) <- c("chrom","start", "end", "rsid", "index_snp")
processed_files <- list.files(pattern="*.processed2$", recursive=F)

my_colnames <- c("table_id","study_id","index_SNP_rsid","current_SNP_rsid","gene_name", "beta", "beta_se",  "p.value","FDR","posterior_prob","Bayes_factor","score","sig_threshold","effect_allele",  "tissue", "sample_size","study_type", "cis/trans","evidence_weight",  "sig_values_no_snp",  "sig_values_no_gene", "adjusted_score")

loadFile <- function(x) {
  print (x)
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  colnames(df) <- my_colnames
  #Remove trailing whitespace
  df$table_id <- str_trim(df$table_id)
  df$study_id <- str_trim(df$study_id)
  df$index_SNP_rsid <- str_trim(df$index_SNP_rsid)
  df$current_SNP_rsid <- str_trim(df$current_SNP_rsid)
  df$gene_name <- str_trim(df$gene_name)
  df$effect_allele <- str_trim(df$effect_allele)
  df$tissue <- str_trim(df$tissue)
  df$study_type <- str_trim(df$study_type)
  df$`cis/trans` <- str_trim(df$`cis/trans`)
  df
}

all_normal <- lapply(processed_files, loadFile)
all_normal_together <- do.call(rbind,all_normal)
all_normal_together <- as.data.frame(all_normal_together)
rm(all_normal)

all_normal_together_tbl <- tbl_df(all_normal_together)
rm(all_normal_together)

#Attach number of experiments per study to the main tibble.
number_of_exp <- all_normal_together_tbl %>% group_by(study_id) %>% summarise(n_of_exps = n_distinct(table_id)) 
all_normal_together_tbl <- left_join(all_normal_together_tbl, number_of_exp, by = "study_id")

#Fix errors in Wang CSRE and Javierre2016_active_promoter_enhancer
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(tissue = replace(tissue, table_id == "Wang2017_csre", study_id[table_id == "Wang2017_csre"]))
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(study_id = replace(study_id, table_id == "Wang2017_csre", 'CSRE'))
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(tissue = replace(tissue, table_id == "Javierre2016_active_promoter_enhancer", study_id[table_id == "Javierre2016_active_promoter_enhancer"]))
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(study_id = replace(study_id, table_id == "Javierre2016_active_promoter_enhancer", 'Javierre2016'))

#Rename categories to more general
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(study_type = recode(study_type, 'oe' = 'promoter-enhancer', 'MTGDR classifier' = 'DGE', 'promoter enhancer'='promoter-enhancer', 'enhancer'='promoter-enhancer','y'='promoter-enhancer','enhancers'='promoter-enhancer','bait'='promoter-enhancer', 'DGE meta'='DGE', 'pir'='promoter-enhancer','x'='promoter-enhancer','regulatory variant predictions'='regulatory variant prediction', 'enhancer-promoter'='promoter-enhancer', 'sQTL'='eQTL', 'coloc - SNP'='coloc', 'exon-ratio eQTL'='eQTL', 'polyA ratio level eQTL'='eQTL', 'FAIRE'='FAIRE-seq', 'isoQTL'='eQTL', 'DEG'='DGE', 'CTCF sites'='CTCF', 'aseQTL'='eQTL', 'coloc - gene'='coloc', 'exon eQTL'='eQTL', 'context specific eQTL'='eQTL', 'ase hQTL'='hQTL', 'transcript ratio eQTL'='eQTL', 'network analysis'='PrixFixe', 'DNMT3B_peaks' = "ChIP-seq", 'ChIP-seq peaks'='ChIP-Seq', 'DNMT3A_peaks'='ChIP-Seq',  'epigenome'='active chromatin state', 'miRNA target sites'='miRNA', 'miRNA seed regions'='miRNA','ChIP-Seq'='ChIP-seq', 'Converved genomic regulatory blocks'='Conserved genomic regulatory blocks', 'ChIP-Seq'='ChIP-seq'))

write.table(all_normal_together_tbl, "all_combined.txt", quote=FALSE, sep="\t", row.names=F)

#all_normal_together_tbl$n_of_exps <- as.numeric(all_normal_together_tbl$n_of_exps)
#all_normal_together_tbl$beta <- as.numeric(all_normal_together_tbl$beta)
#all_normal_together_tbl$beta_se <- as.numeric(all_normal_together_tbl$beta_se)
#all_normal_together_tbl$p.value <- as.numeric(all_normal_together_tbl$p.value)
#all_normal_together_tbl$FDR <- as.numeric(all_normal_together_tbl$FDR)
#all_normal_together_tbl$posterior_prob <- as.numeric(all_normal_together_tbl$posterior_prob)
#all_normal_together_tbl$Bayes_factor <- as.numeric(all_normal_together_tbl$Bayes_factor)
#all_normal_together_tbl$score <- as.numeric(all_normal_together_tbl$score)
#all_normal_together_tbl$evidence_weight <- as.numeric(all_normal_together_tbl$evidence_weight)
#all_normal_together_tbl$sig_values_no_snp <- as.numeric(all_normal_together_tbl$sig_values_no_snp)
#all_normal_together_tbl$sig_values_no_gene <- as.numeric(all_normal_together_tbl$sig_values_no_gene)
#all_normal_together_tbl$adjusted_score <- as.numeric(all_normal_together_tbl$adjusted_score)
#all_normal_together_tbl$sig_threshold <- as.numeric(all_normal_together_tbl$sig_threshold)
#all_normal_together_tbl$sample_size <- as.numeric(all_normal_together_tbl$sample_size)
#write.table(all_normal_together_tbl, "all_processed2.tbl", quote=F, sep="\t", row.names=F)

#all_normal_together <- read.delim("all_processed2.tbl", header=T, stringsAsFactors=T,row#.names=NULL, sep="\t")



#Model 00 Rank genes by counts, divide by index locus and print. This prioritizes genes with matches to many SNPs. 
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(adjusted_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model00_gene.txt", quote=F, sep="\t")
rank_snps <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(adjusted_score))
names(rank_snps)[3] <- "final_score"
rank_snps <- rank_snps %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model00_snp.txt", quote=F, sep="\t")
#Model_05
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=adjusted_score/sqrt(n_of_exps))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model05_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=adjusted_score/sqrt(n_of_exps))
rank_snps <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps)[3] <- "final_score"
rank_snps <- rank_snps %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model05_snp.txt", quote=F, sep="\t")
#Model_06
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=(adjusted_score/sqrt(n_of_exps)/sqrt(sig_values_no_gene)))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model06_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=(adjusted_score/sqrt(n_of_exps)/sqrt(sig_values_no_snp)))
rank_snps <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps)[3] <- "final_score"
rank_snps <- rank_snps %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model06_snp.txt", quote=F, sep="\t")
#Model_07
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, adjusted_score/sqrt(n_of_exps), ifelse(evidence_weight==2, 2*adjusted_score/sqrt(n_of_exps), 10*adjusted_score/sqrt(n_of_exps))))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model07_gene.txt", quote=F, sep="\t")

rank_snps <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps)[3] <- "final_score"
rank_snps <- rank_snps %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model07_snp.txt", quote=F, sep="\t")
#Model_08
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps))/sig_values_no_gene, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene, 10*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene)))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model08_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, 10*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp)))
rank_snps <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps)[3] <- "final_score"
rank_snps <- rank_snps %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model08_snp.txt", quote=F, sep="\t")
#Model_09
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, adjusted_score/sqrt(n_of_exps), ifelse(evidence_weight==2, 2*adjusted_score/sqrt(n_of_exps), 20*adjusted_score/sqrt(n_of_exps))))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model09_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, adjusted_score/sqrt(n_of_exps), ifelse(evidence_weight==2, 2*adjusted_score/sqrt(n_of_exps), 20*adjusted_score/sqrt(n_of_exps))))
rank_snps_09 <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps_09)[3] <- "final_score"
rank_snps_09 <- rank_snps_09 %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps_09, "Model09_snp.txt", quote=F, sep="\t")
#Model_10
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps))/sig_values_no_gene, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene, 20*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene)))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model10_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, 20*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp)))
rank_snps_10 <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(final_score))
names(rank_snps_10)[3] <- "final_score"
rank_snps_10 <- rank_snps_10 %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps_10, "Model10_snp.txt", quote=F, sep="\t")

#Model11
#Adjust for square root of number of SNPs
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene / sqrt(sig_values_no_snp), ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene / sqrt(sig_values_no_snp), 20*(adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene / sqrt(sig_values_no_snp))))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model11_gene.txt", quote=F, sep="\t")

#Model12
#Adjust for square root of number of SNPs and increase weight for evidence 1
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene / sqrt(sig_values_no_snp), ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene / sqrt(sig_values_no_snp), 30*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene / sqrt(sig_values_no_snp))))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model12_gene.txt", quote=F, sep="\t")

#Model13
#Adjust for square root of number of SNPs and prioritize genes wih associated with lots of study types/study ids
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, ( adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene / sqrt(sig_values_no_snp), ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene / sqrt(sig_values_no_snp), 20*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene / sqrt(sig_values_no_snp))))
rank_genes_13 <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(sqrt((length(unique(study_type)) + length(unique(study_id)))/2)*final_score))
names(rank_genes_13)[3] <- "final_score"
rank_genes_13 <- rank_genes_13 %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes_13, "Model13_gene.txt", quote=F, sep="\t")

#Model14
#Prioritize genes wih associated with lots of study types/study ids
all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps)) / sig_values_no_gene, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene, 20*(adjusted_score/sqrt(n_of_exps))/sig_values_no_gene)))
rank_genes <- all_normal_together_tbl %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sum(sqrt((length(unique(study_type)) + length(unique(study_id)))/2)*final_score))
names(rank_genes)[3] <- "final_score"
rank_genes <- rank_genes %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% filter(!is.na(final_score))
write.table(rank_genes, "Model14_gene.txt", quote=F, sep="\t")

all_normal_together_tbl <- all_normal_together_tbl %>% mutate(final_score=ifelse(evidence_weight==3, (adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, ifelse(evidence_weight==2, 2*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp, 20*(adjusted_score/sqrt(n_of_exps))/sig_values_no_snp)))
rank_snps_14 <- all_normal_together_tbl %>% group_by(index_SNP_rsid, current_SNP_rsid) %>% summarise(sum(sqrt((length(unique(study_type)) + length(unique(study_id)))/2)*final_score))
names(rank_snps_14)[3] <- "final_score"
rank_snps_14 <- rank_snps_14 %>% arrange(desc(final_score)) %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% filter(!is.na(final_score))
write.table(rank_snps, "Model14_snp.txt", quote=F, sep="\t")


#5 Print all evidence (matching rows) for a given gene for a given locus
for (f in unique(all_normal_together_tbl$index_SNP_rsid)) {
  my_output_0 <- all_normal_together_tbl %>% filter(index_SNP_rsid == f)
  for (z in my_output_0$gene_name)
  {
    my_output <- my_output_0 %>% filter(gene_name== z)
    my_filename <- paste(f, "_", z, sep="", ".genes")
    write.table(my_output, my_filename, quote=F, sep="\t", row.names=F)  
  }}
#5 Print all evidence (matching rows) for a given SNP for a given locus
for (f in unique(all_normal_together$index_SNP_rsid)) {
  my_output_0 <- all_normal_together_tbl %>% filter(index_SNP_rsid == f)
  for (z in my_output_0$current_SNP_rsid)
  {
    my_output <- my_output_0 %>% filter(current_SNP_rsid == z)
    my_filename <- paste(f, "_", z, sep="", ".snps")
    write.table(my_output, my_filename, quote=F, sep="\t", row.names=F)   
  }}

#Print all evidence for top 10 genes for Model 13 (genes)
for (f in unique(all_normal_together_tbl$index_SNP_rsid)) {
ranked <- rank_genes_13 %>% filter(index_SNP_rsid == f)
top10_genes <- head(ranked, 10)
for (z in top10_genes$gene_name)
{
  my_output <- all_normal_together_tbl %>% filter(index_SNP_rsid == f & gene_name== z)
  my_filename <- paste(f, "_", z, sep="", ".genes")
  my_path <- paste("./Model13_top_10_gene/", my_filename, sep="")
  write.table(my_output, my_path, quote=F, sep="\t", row.names=F)  
}}

#Print all evidence for top 10 SNPs for Model 10 (SNPs)
for (f in unique(all_normal_together_tbl$index_SNP_rsid)) {
  ranked <- rank_snps_10 %>% filter(index_SNP_rsid == f)
  top10_snps <- head(ranked, 10)
  for (z in top10_snps$current_SNP_rsid)
  {
    my_output <- all_normal_together_tbl %>% filter(index_SNP_rsid == f & current_SNP_rsid == z)
    my_filename <- paste(f, "_", z, sep="", ".snps")
    my_path <- paste("./Model10_top_10_snp/", my_filename, sep="")
    write.table(my_output, my_path, quote=F, sep="\t", row.names=F)   
  }}

#Print all evidence for top 10 SNPs for Model 9 (SNPs)
for (f in unique(all_normal_together_tbl$index_SNP_rsid)) {
  ranked <- rank_snps_09 %>% filter(index_SNP_rsid == f)
  top10_snps <- head(ranked, 10)
  for (z in top10_snps$current_SNP_rsid)
  {
    my_output <- all_normal_together_tbl %>% filter(index_SNP_rsid == f & current_SNP_rsid == z)
    my_filename <- paste(f, "_", z, sep="", ".snps")
    my_path <- paste("./Model09_top_10_snp/", my_filename, sep="")
    write.table(my_output, my_path, quote=F, sep="\t", row.names=F)   
  }}

#Print all evidence for top 10 SNPs for Model 14 (SNPs)
for (f in unique(all_normal_together_tbl$index_SNP_rsid)) {
  ranked <- rank_snps_14 %>% filter(index_SNP_rsid == f)
  top10_snps <- head(ranked, 10)
  for (z in top10_snps$current_SNP_rsid)
  {
    my_output <- all_normal_together_tbl %>% filter(index_SNP_rsid == f & current_SNP_rsid == z)
    my_filename <- paste(f, "_", z, sep="", ".snps")
    my_path <- paste("./Model14_top_10_snp/", my_filename, sep="")
    write.table(my_output, my_path, quote=F, sep="\t", row.names=F)   
  }}

#12 Separate genes which are in the 3Mbp list and the genes beyond them in #5 for ranking and printing. Lookup the second category - annotate their start and end, see how far they are from index SNP and print a second table, similar to paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.hugo.
genes_in_df <- levels(all_normal_together_tbl$gene_name)
genes_in_3mbp_interval <- my_gene_ref$ensembl
only_present_in_df <- setdiff(genes_in_df,genes_in_3mbp_interval)
only_present_in_interval <- setdiff(genes_in_3mbp_interval,genes_in_df)
library("rtracklayer")
my_human <- import.gff("../Homo_sapiens.GRCh37.87.gff3.gz")
engs <- grep("^ENSG", only_present_in_df)
my_engs_names_to_lookup <- only_present_in_df[engs]
my_reg_names_to_lookup <- toupper(only_present_in_df[-engs])

my_eng_matches <- match(my_engs_names_to_lookup, my_human$gene_id)
my_human$UppercaseName <- toupper(my_human$Name)
my_reg_matches <- match(my_reg_names_to_lookup,my_human$UppercaseName)

my_eng_ranges <- my_human[my_eng_matches[!is.na(my_eng_matches)]]
my_reg_ranges <- my_human[my_reg_matches[!is.na(my_reg_matches)]]

df <- data.frame(chr=seqnames(my_reg_ranges),
                 starts=start(my_reg_ranges),
                 ends=end(my_reg_ranges),
                 ensembl=my_reg_ranges$gene_id,
                 gene_name=my_reg_ranges$Name)
df2 <- data.frame(chr=seqnames(my_eng_ranges),
                  starts=start(my_eng_ranges),
                  ends=end(my_eng_ranges),
                  ensembl=my_eng_ranges$gene_id,
                  gene_name=my_eng_ranges$Name)

secondary_alltogether <- rbind(df, df2)
secondary_alltogether <- tbl_df(secondary_alltogether)

#Add the index SNP to each gene to know which locus we are focussing on.
match_index <- match(secondary_alltogether$gene_name, all_normal_together_tbl$gene_name)
match_index_subset <- all_normal_together_tbl[match_index,]
secondary_alltogether$snp_name <- match_index_subset$index_SNP_rsid
secondary_alltogether <- secondary_alltogether %>% select("chr", "starts", "ends", "snp_name", "ensembl", "gene_name") %>% arrange(chr,starts)

write.table(secondary_alltogether, "paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches.secondary", quote=FALSE, row.names=F, sep="\t")

