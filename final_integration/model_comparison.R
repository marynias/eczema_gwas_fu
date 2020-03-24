list_index_SNPs <- read.delim("index_SNPs.txt", header=T, sep="\t", stringsAsFactors=F)
closest_genes <- read.delim("closest_genes.txt", header=T, sep="\t", stringsAsFactors=F)
best_gene <-read.delim("best_genes.txt", header=T, sep="\t", stringsAsFactors=F)
best_snp <-read.delim("best_SNPs.txt", header=T, sep="\t", stringsAsFactors=F)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
best_gene$Gene <- str_trim(best_gene$Gene)
#Read in a list of models.
model_genes <- list.files(pattern="*_gene_model.txt$", recursive=F)
model_snps <- list.files(pattern="*_snp_model.txt$", recursive=F)
loadFile <- function(x) {
  print (x)
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  df <- tbl_df(df)
}

add_rank <- function(x) {
x %>% group_by(index_SNP_rsid) %>% mutate(rank = min_rank(desc(final_score))) %>%
                                            ungroup()
}

extract_top <- function(x, my_n) {
module <- split(x, x$index_SNP_rsid)
df <- lapply(module, head, n=my_n)
}

all_gene_models <- lapply(model_genes, loadFile)
all_snp_models <- lapply(model_snps, loadFile)

all_gene_models <- lapply(all_gene_models, add_rank)
all_snp_models <- lapply(all_snp_models, add_rank)

#Print ranked tables to file
for(i in seq_along(all_gene_models)){
out_name = paste(all_gene_models[[i]]$Model[1], "_gene_ranked.txt", sep="")
write.table(all_gene_models[[i]][c("index_SNP_rsid", "gene_name", "final_score", "rank")], out_name, quote=F, sep="\t", row.names=F)
}

#Print ranked tables to file
for(i in seq_along(all_snp_models)){
  out_name = paste(all_snp_models[[i]]$Model[1], "_snp_ranked.txt", sep="")
  write.table(all_snp_models[[i]][c("index_SNP_rsid", "current_SNP_rsid", "final_score", "rank")], out_name, quote=F, sep="\t", row.names=F)
}


#Look at presence of index loci among the SNP ranks for each locus.
res_snp <- vector('list', length(all_snp_models))
for(i in seq_along(all_snp_models)){
  res_snp[[i]] <- all_snp_models[[i]][all_snp_models[[i]]$index_SNP_rsid==all_snp_models[[i]]$current_SNP_rsid,]
}

res_snp_select <- vector('list', length(all_snp_models))
for(i in seq_along(all_snp_models)){
  res_snp_select[[i]] <- all_snp_models[[i]][all_snp_models[[i]]$index_SNP_rsid %in% best_snp$Top_SNP & all_snp_models[[i]]$current_SNP_rsid %in% best_snp$Top_SNP,]
}

res_snp_names <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model14")
#Filter gene tibble to keep only genes with rank smaller than 6.
for(i in seq_along(res_snp)){
  res_snp[[i]] <- res_snp[[i]] %>% filter(rank < 6) %>% mutate(Model=res_snp_names[i])
}
for(i in seq_along(res_snp_select)){
  res_snp_select[[i]] <- res_snp_select[[i]] %>% filter(rank < 6) %>% mutate(Model=res_snp_names[i])
}


combined_df2_snp <- do.call(rbind, res_snp)
combined_df2_snp_select <- do.call(rbind, res_snp_select)
#Write the SNP rank to file.
combined_df2_snp <- combined_df2_snp[c("index_SNP_rsid", "current_SNP_rsid", "final_score", "rank", "Model")]
combined_df2_snp_select <- combined_df2_snp_select[c("index_SNP_rsid", "current_SNP_rsid", "final_score", "rank", "Model")]
write.table(combined_df2_snp, "snp_ranking.txt", quote=F, sep="\t", row.names=F)
write.table(combined_df2_snp_select, "snp_ranking_select.txt", quote=F, sep="\t", row.names=F)

snp_ranking_tbl = tbl_df(combined_df2_snp)
snp_ranking_tbl_select = tbl_df(combined_df2_snp_select)

#All the loci
snp_ranking_summary_all <- snp_ranking_tbl %>% group_by(Model) %>% summarize(total = n(), average_rank = mean(rank)) %>% arrange(desc(total))
write.table(snp_ranking_summary_all, "all_snps_model_rank.txt", quote=F, sep="\t", row.names=F)

#Only the top hits.
snp_ranking_summary_select <- snp_ranking_tbl_select %>% group_by(Model) %>% summarize(total = n(), average_rank = mean(rank)) %>% arrange(desc(total))
write.table(snp_ranking_summary_select, "select_snps_model_rank.txt", quote=F, sep="\t", row.names=F)

#Convert to long format for plotting
melted_snp <- as.data.frame(melt(snp_ranking_summary_all, id.vars=c("Model")))
melted_snp2 <- as.data.frame(melt(snp_ranking_summary_select, id.vars=c("Model")))

theme_set(theme_gray(base_size = 16))
ggplot() + geom_col(data=melted_snp, aes(x = Model, y = value, fill = variable), position = "dodge") + theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_manual(values=c("black", "red2")) + ggtitle ("Presence of GWAS index SNPs among top 5 SNPs per locus")

theme_set(theme_gray(base_size = 16))
ggplot() + geom_col(data=melted_snp2, aes(x = Model, y = value, fill = variable), position = "dodge") + theme(legend.title=element_blank()) + scale_fill_manual(values=c("black", "red2")) + ggtitle ("Presence of GWAS index SNPs among top 5 SNPs in 5 select loci")

#Look at presence of nearby genes among the gene ranks for each locus.
res_gene <- vector('list', length(all_gene_models))
for(i in seq_along(all_gene_models)){
res_gene[[i]] <- all_gene_models[[i]][all_gene_models[[i]]$index_SNP_rsid %in% closest_genes$Top_SNP & all_gene_models[[i]]$gene_name %in% closest_genes$Gene,]
}

res_gene_select <- vector('list', length(all_gene_models))
for(i in seq_along(all_gene_models)){
  res_gene_select[[i]] <- all_gene_models[[i]][all_gene_models[[i]]$index_SNP_rsid %in% best_gene$Top_SNP & all_gene_models[[i]]$gene_name %in% best_gene$Gene,]
}


res_gene_names <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model11", "Model12", "Model13", "Model14")
#Filter gene tibble to keep only genes with rank smaller than 6.
for(i in seq_along(res_gene)){
  res_gene[[i]] <- res_gene[[i]] %>% filter(rank < 6) %>% mutate(Model=res_gene_names[i])
}
for(i in seq_along(res_gene_select)){
  res_gene_select[[i]] <- res_gene_select[[i]] %>% filter(rank < 6) %>% mutate(Model=res_gene_names[i])
}

combined_df2 <- do.call(rbind, res_gene)
combined_df2_select <- do.call(rbind, res_gene_select)
#Write the gene rank to file.
combined_df2 <- combined_df2[c("index_SNP_rsid", "gene_name", "final_score", "rank", "Model")]
combined_df2_select <- combined_df2_select[c("index_SNP_rsid", "gene_name", "final_score", "rank", "Model")]
write.table(combined_df2, "gene_ranking.txt", quote=F, sep="\t", row.names=F)
write.table(combined_df2_select, "gene_ranking_select.txt", quote=F, sep="\t", row.names=F)

gene_ranking_tbl = tbl_df(combined_df2)
gene_ranking_tbl_select = tbl_df(combined_df2_select)


#All the loci
gene_ranking_summary_all <- gene_ranking_tbl %>% group_by(Model) %>% summarize(total = n(), average_rank = mean(rank)) %>% arrange(desc(total))
write.table(gene_ranking_summary_all, "all_genes_model_rank.txt", quote=F, sep="\t", row.names=F)

#Only the top hits.
gene_ranking_summary_select <- gene_ranking_tbl_select %>% group_by(Model) %>% summarize(total = n(), average_rank = mean(rank)) %>% arrange(desc(total))
write.table(gene_ranking_summary_select, "select_genes_model_rank.txt", quote=F, sep="\t", row.names=F)

#Convert to long format for plotting
melted_gene <- as.data.frame(melt(gene_ranking_summary_all, id.vars=c("Model")))
melted_gene2 <- as.data.frame(melt(gene_ranking_summary_select, id.vars=c("Model")))

theme_set(theme_gray(base_size = 16))
ggplot() + geom_col(data=melted_gene, aes(x = Model, y = value, fill = variable), position = "dodge") + theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_manual(values=c("black", "blue")) + ggtitle ("Presence of closest genes among top 5 genes per locus")

theme_set(theme_gray(base_size = 16))
ggplot() + geom_col(data=melted_gene2, aes(x = Model, y = value, fill = variable), position = "dodge") + theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_manual(values=c("black", "blue")) + ggtitle ("Presence of best genes among top 5 genes per locus")

#Plot overall score for all intervals
ggplot(toy_model_gene_tbl, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(toy_model_gene_tbl, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,2000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model10_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(toy_model_snp_tbl, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(toy_model_snp_tbl, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,50)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model10_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(toy_model_gene_tbl2, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(toy_model_gene_tbl2, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,30000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model09_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(toy_model_snp_tbl2, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(toy_model_snp_tbl2, aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,700)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model09_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[1]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[1]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,30000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model00_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[1]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[1]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,600)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model00_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[2]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[2]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,13000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model05_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[2]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[2]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,400)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model05_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[3]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[3]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,5000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model06_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[3]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[3]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,20)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model06_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[4]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[4]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,30000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model07_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[4]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[4]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,700)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model07_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[5]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[5]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,2000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) ggsave("Model08_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[5]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[5]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,25)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model08_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_snp_models[[8]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_snp_models[[8]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,100)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model14_snp.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[8]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[8]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,200)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model11_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[9]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[9]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,200)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model12_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[10]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[10]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,500)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model13_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

ggplot(all_gene_models[[11]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA))
ggplot(all_gene_models[[11]], aes(x=final_score)) + geom_histogram(bins=100, fill="aquamarine3") + coord_cartesian(xlim=c(0,1000)) + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave("Model14_gene.pdf", plot=last_plot(), width = 80, height = 8, limitsize = FALSE)

#Plot top 5 score for all intervals.
plot_top_snp_model <- function(x,y) {
  my_input <- x[x$rank <6,] %>% arrange(desc(final_score))
my_plot <-  ggplot(my_input, aes(y=final_score, x=current_SNP_rsid)) +  geom_bar(stat="identity",fill="palevioletred") + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
ggsave(y, plot=last_plot(), width = 80, height = 8, limitsize = FALSE)  
}

plot_top_snp_model(all_snp_models[[1]],"Model00_top_snp.pdf")
plot_top_snp_model(all_snp_models[[2]],"Model05_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[3]],"Model06_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[4]],"Model07_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[5]],"Model08_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[6]],"Model09_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[7]],"Model10_top_snp.pdf" )
plot_top_snp_model(all_snp_models[[8]],"Model14_top_snp.pdf" )

plot_top_gene_model <- function(x,y) {
  my_input <- x[x$rank <6,] %>% arrange(desc(final_score))
  my_plot <-  ggplot(my_input, aes(y=final_score, x=gene_name)) +  geom_bar(stat="identity",fill="palevioletred") + facet_grid(. ~ index_SNP_rsid) + theme(panel.grid.major.y = element_line(colour = "black",size=0.2), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA)) 
  ggsave(y, plot=last_plot(), width = 80, height = 8, limitsize = FALSE)  
}
plot_top_gene_model(all_gene_models[[1]],"Model00_top_gene.pdf")
plot_top_gene_model(all_gene_models[[2]],"Model05_top_gene.pdf")
plot_top_gene_model(all_gene_models[[3]],"Model06_top_gene.pdf")
plot_top_gene_model(all_gene_models[[4]],"Model07_top_gene.pdf")
plot_top_gene_model(all_gene_models[[5]],"Model08_top_gene.pdf")
plot_top_gene_model(all_gene_models[[6]],"Model09_top_gene.pdf")
plot_top_gene_model(all_gene_models[[7]],"Model10_top_gene.pdf")
plot_top_gene_model(all_gene_models[[8]],"Model11_top_gene.pdf")
plot_top_gene_model(all_gene_models[[9]],"Model12_top_gene.pdf")
plot_top_gene_model(all_gene_models[[10]],"Model13_top_gene.pdf")
plot_top_gene_model(all_gene_models[[11]],"Model14_top_gene.pdf")

correlation_top_snp <- data.frame(matrix(nrow=length(all_snp_models), ncol=length(all_snp_models))) 
sd_top_snp <- data.frame(matrix(nrow=length(all_snp_models), ncol=length(all_snp_models))) 
for (i in seq_along(all_snp_models)){
  for (j in seq_along(all_snp_models)) {
    rho_values = double()
    for (my_index in unique(all_snp_models[[i]]$index_SNP_rsid)) {
      #Sort by rsid
      top_snp1 <- all_snp_models[[i]][all_snp_models[[i]]$index_SNP_rsid==my_index & all_snp_models[[i]]$rank <11,]
      top_snp2 <- all_snp_models[[j]][all_snp_models[[j]]$index_SNP_rsid==my_index & all_snp_models[[j]]$rank <11,] 
      all_top <- unique(c(top_snp1$current_SNP_rsid, top_snp2$current_SNP_rsid))
      sorted1 <- all_snp_models[[i]][all_snp_models[[i]]$index_SNP_rsid==my_index & all_snp_models[[i]]$current_SNP_rsid %in% all_top,] %>% arrange(current_SNP_rsid)
      sorted2 <- all_snp_models[[j]][all_snp_models[[j]]$index_SNP_rsid==my_index & all_snp_models[[j]]$current_SNP_rsid %in% all_top,] %>% arrange(current_SNP_rsid)
      out <- cor.test(sorted1$final_score, sorted2$final_score, method = "spearman", exact=FALSE) 
      my_rho <-  out$estimate
      rho_values <- c(rho_values, my_rho)
    }
    mean_rho <- mean(rho_values)
    sd_rho <- sd(rho_values)
    correlation_top_snp[i, j] <- mean_rho
    sd_top_snp[i, j] <- sd_rho
    
  }
}


correlation_all_snp <- data.frame(matrix(nrow=length(all_snp_models), ncol=length(all_snp_models))) 
sd_all_snp <- data.frame(matrix(nrow=length(all_snp_models), ncol=length(all_snp_models))) 
for (i in seq_along(all_snp_models)){
  for (j in seq_along(all_snp_models)) {
    rho_values = double()
    for (my_index in unique(all_snp_models[[i]]$index_SNP_rsid)) {
      #Sort by rsid
      sorted1 <- all_snp_models[[i]][all_snp_models[[i]]$index_SNP_rsid==my_index,] %>% arrange(current_SNP_rsid)
      sorted2 <- all_snp_models[[j]][all_snp_models[[j]]$index_SNP_rsid==my_index,] %>% arrange(current_SNP_rsid)
      out <- cor.test(sorted1$final_score, sorted2$final_score, method = "spearman", exact=FALSE) 
      my_rho <-  out$estimate
      rho_values <- c(rho_values, my_rho)
    }
    mean_rho <- mean(rho_values)
    sd_rho <- sd(rho_values)
    correlation_all_snp[i, j] <- mean_rho
    sd_all_snp[i, j] <- sd_rho
    
  }
}

correlation_top_gene <- data.frame(matrix(nrow=length(all_gene_models), ncol=length(all_gene_models)))
sd_top_gene <- data.frame(matrix(nrow=length(all_gene_models), ncol=length(all_gene_models))) 
for (i in seq_along(all_gene_models)){
  for (j in seq_along(all_gene_models)) {
    rho_values = double()
    for (my_index in unique(all_gene_models[[i]]$index_SNP_rsid)) {
      #Sort by rsid
      top_gene1 <- all_gene_models[[i]][all_gene_models[[i]]$index_SNP_rsid==my_index & all_gene_models[[i]]$rank <11,]
      top_gene2 <- all_gene_models[[j]][all_gene_models[[j]]$index_SNP_rsid==my_index & all_gene_models[[j]]$rank <11,] 
      all_top <- unique(c(top_gene1$gene_name, top_gene2$gene_name))
      sorted1 <- all_gene_models[[i]][all_gene_models[[i]]$index_SNP_rsid==my_index & all_gene_models[[i]]$gene_name %in% all_top,] %>% arrange(gene_name)
      sorted2 <- all_gene_models[[j]][all_gene_models[[j]]$index_SNP_rsid==my_index & all_gene_models[[j]]$gene_name %in% all_top,] %>% arrange(gene_name)
      out <- cor.test(sorted1$final_score, sorted2$final_score, method = "spearman", exact=FALSE) 
      my_rho <-  out$estimate
      rho_values <- c(rho_values, my_rho)
    }
    mean_rho <- mean(rho_values)
    sd_rho <- sd(rho_values)
    correlation_top_gene[i, j] <- mean_rho
    sd_top_gene[i, j] <- sd_rho
    
  }
}

correlation_all_gene <- data.frame(matrix(nrow=length(all_gene_models), ncol=length(all_gene_models)))
sd_all_gene <- data.frame(matrix(nrow=length(all_gene_models), ncol=length(all_gene_models))) 
for (i in seq_along(all_gene_models)){
  for (j in seq_along(all_gene_models)) {
    rho_values = double()
    for (my_index in unique(all_gene_models[[i]]$index_SNP_rsid)) {
      #Sort by rsid
      sorted1 <- all_gene_models[[i]][all_gene_models[[i]]$index_SNP_rsid==my_index,] %>% arrange(gene_name)
      sorted2 <- all_gene_models[[j]][all_gene_models[[j]]$index_SNP_rsid==my_index,] %>% arrange(gene_name)
      out <- cor.test(sorted1$final_score, sorted2$final_score, method = "spearman", exact=FALSE) 
      my_rho <-  out$estimate
      rho_values <- c(rho_values, my_rho)
    }
    mean_rho <- mean(rho_values)
    sd_rho <- sd(rho_values)
    correlation_all_gene[i, j] <- mean_rho
    sd_all_gene[i, j] <- sd_rho
    
  }
}



format_and_write_out <- function(x, y) {
  colnames(x) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model11", "Model12", "Model13", "Model14")
  rownames(x) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model11", "Model12", "Model13", "Model14")
  write.table(format(x,digits=2), y, sep='\t',row.names=T, quote=F, col.names = NA)  
  x  
}

format_and_write_out2 <- function(x, y) {
  colnames(x) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model14")
  rownames(x) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model14")
  write.table(format(x,digits=2), y, sep='\t',row.names=T, quote=F, col.names = NA)  
  x  
}

correlation_top_snp <- format_and_write_out2(correlation_top_snp, "correlation_top_snp.txt")
correlation_all_snp <- format_and_write_out2(correlation_all_snp, "correlation_all_snp.txt")
correlation_top_gene <- format_and_write_out(correlation_top_gene, "correlation_top_gene.txt")
correlation_all_gene <- format_and_write_out(correlation_all_gene, "correlation_all_gene.txt")
sd_top_snp <- format_and_write_out2(sd_top_snp, "sd_top_snp.txt")
sd_all_snp <- format_and_write_out2(sd_all_snp, "sd_all_snp.txt")
sd_top_gene <- format_and_write_out(sd_top_gene, "sd_top_gene.txt")
sd_all_gene <- format_and_write_out(sd_all_gene, "sd_all_gene.txt")

plot_coloured <- function(x) {
  x$Model <- rownames(x)
  plotDat <- gather(x, key = "Model2", value = "rho", -Model)
  plotDat$rho <- round(plotDat$rho,2)
  
  ggplot(plotDat, aes(Model2, Model, col = rho, fill = rho, label = rho)) +
    geom_tile() +
    geom_text(col = "black") +
    theme_minimal() + theme(axis.title=element_blank()) +
    scale_fill_gradient2(midpoint=0,low = "blue", mid = "white", high = "red") +
    scale_color_gradient2(midpoint=0,low = "blue", mid = "white", high = "red")  
  
}

plot_coloured2 <- function(x) {
  x$Model <- rownames(x)
  plotDat <- gather(x, key = "Model2", value = "rho", -Model)
  plotDat$rho <- round(plotDat$rho,2)
  
  ggplot(plotDat, aes(Model2, Model, col = rho, fill = rho, label = rho)) +
    geom_tile() +
    geom_text(col = "black") +
    theme_minimal() + theme(axis.title=element_blank()) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red", limits=c(0, 1)) +
    scale_color_gradient2(low = "white", mid = "yellow", high = "red", limits=c(0, 1))  
  
}


plot_coloured(correlation_top_snp)
plot_coloured(correlation_all_snp)
plot_coloured2(correlation_top_gene)
plot_coloured2(correlation_all_gene)
plot_coloured2(sd_top_gene)
plot_coloured2(sd_all_gene)
plot_coloured2(sd_all_snp)
plot_coloured2(sd_top_snp)

names(all_gene_models) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model11", "Model12", "Model13", "Model14")
all_gene_models[[1]]$Model <- "Model00"
all_gene_models[[2]]$Model <- "Model05"
all_gene_models[[3]]$Model <- "Model06"
all_gene_models[[4]]$Model <- "Model07"
all_gene_models[[5]]$Model <- "Model08"
all_gene_models[[6]]$Model <- "Model09"
all_gene_models[[7]]$Model <- "Model10"
all_gene_models[[8]]$Model <- "Model11"
all_gene_models[[9]]$Model <- "Model12"
all_gene_models[[10]]$Model <- "Model13"
all_gene_models[[11]]$Model <- "Model14"

pick_genes <- function(x, y,z) {
sub_df <- x[x$index_SNP_rsid==my_index & x$gene_name %in% out,]  
transposed_mini<- t(sub_df)
colnames(transposed_mini) <- transposed_mini[3,]
transposed_mini <- as.data.frame(transposed_mini)
transposed_mini$Model <- transposed_mini[6,1]
final_mini <- transposed_mini[5,]
rownames(final_mini) <- final_mini$Model
final_mini$Model <- NULL
final_mini
}

plot_coloured_rank_gene <- function(x,y) {
  x$Model <- rownames(x)
  plotDat <- gather(x, key = "Model2", value = "rank", -Model)
  plotDat$inverse <- 1/plotDat$rank
  
  ggplot(plotDat, aes(Model2, Model, col = inverse, fill = inverse, label = rank)) +
    geom_tile() +
    geom_text(col = "black") +
    theme_minimal() + theme(axis.title=element_blank()) +
    scale_fill_gradient2(low = "white", mid="gray", midpoint=1/8, high = "green") + scale_color_gradient2(low = "white", mid="gray", midpoint=1/8, high = "green")
  ggsave(y, plot=last_plot(), width = 7, height = 5, limitsize = FALSE)
  
}


for (my_index in unique(all_gene_models[[1]]$index_SNP_rsid)) {
  out = list()
  for (i in seq_along(all_gene_models)) {
  my_top_5 <- all_gene_models[[i]][all_gene_models[[i]]$rank <6 & all_gene_models[[i]]$index_SNP_rsid==my_index,]$gene_name
  out <- unlist(unique(c(out,  my_top_5)))
  }
  #Select all the genes in the 7 models.
  for_heatmap <- lapply(all_gene_models, pick_genes, out, my_index)
  combined_df_heat <- do.call(rbind, for_heatmap)
  my_out <- paste(my_index, "_gene_top5_all.txt", sep="")
  write.table(t(combined_df_heat), my_out, sep="\t", col.names=NA, quote=F)
  combined_df_heat <- read.delim(my_out, sep="\t", header=TRUE, row.names=1)
  my_out2 <- paste(my_index, "_gene_top5_all.pdf", sep="")
  plot_coloured_rank_gene(combined_df_heat, my_out2)
}

names(all_snp_models) <- c("Model00", "Model05", "Model06", "Model07", "Model08", "Model09", "Model10", "Model14")
all_snp_models[[1]]$Model <- "Model00"
all_snp_models[[2]]$Model <- "Model05"
all_snp_models[[3]]$Model <- "Model06"
all_snp_models[[4]]$Model <- "Model07"
all_snp_models[[5]]$Model <- "Model08"
all_snp_models[[6]]$Model <- "Model09"
all_snp_models[[7]]$Model <- "Model10"
all_snp_models[[8]]$Model <- "Model14"

pick_snps <- function(x, y,z) {
  sub_df <- x[x$index_SNP_rsid==my_index & x$current_SNP_rsid %in% out,]  
  transposed_mini<- t(sub_df)
  colnames(transposed_mini) <- transposed_mini[3,]
  transposed_mini <- as.data.frame(transposed_mini)
  transposed_mini$Model <- transposed_mini[6,1]
  final_mini <- transposed_mini[5,]
  rownames(final_mini) <- final_mini$Model
  final_mini$Model <- NULL
  final_mini
}

plot_coloured_rank_snp <- function(x,y) {
  x$Model <- rownames(x)
  plotDat <- gather(x, key = "Model2", value = "rank", -Model)
  plotDat$inverse <- 1/plotDat$rank
  
  ggplot(plotDat, aes(Model2, Model, col = inverse, fill = inverse, label = rank)) +
    geom_tile() +
    geom_text(col = "black") +
    theme_minimal() + theme(axis.title=element_blank()) +
    scale_fill_gradient2(low = "white", mid="gray", midpoint=1/8, high = "green") + scale_color_gradient2(low = "white", mid="gray", midpoint=1/8, high = "green")
  ggsave(y, plot=last_plot(), width = 7, height = 10, limitsize = FALSE)
  
}

for (my_index in unique(all_snp_models[[1]]$index_SNP_rsid)) {
  out = list()
  for (i in seq_along(all_snp_models)) {
    my_top_5 <- all_snp_models[[i]][all_snp_models[[i]]$rank <6 & all_snp_models[[i]]$index_SNP_rsid==my_index,]$current_SNP_rsid
    out <- unlist(unique(c(out,  my_top_5)))
  }
  #Select all the SNPs in the 7 models.
  for_heatmap <- lapply(all_snp_models, pick_snps, out, my_index)
  combined_df_heat <- do.call(rbind, for_heatmap)
  my_out <- paste(my_index, "_snp_top5_all.txt", sep="")
  write.table(t(combined_df_heat), my_out, sep="\t", col.names=NA, quote=F)
  combined_df_heat <- read.delim(my_out, sep="\t", header=TRUE, row.names=1)
  my_out2 <- paste(my_index, "_snp_top5_all.pdf", sep="")
  plot_coloured_rank_snp(combined_df_heat, my_out2)
}


