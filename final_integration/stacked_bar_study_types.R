library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
#Set my palette

# to get the colors from a palette:
palette1 <- brewer.pal(9,"Set1")
palette1

palette2 <- brewer.pal(8,"Set2")
palette2

palette3 <- brewer.pal(9,"Set3")
palette3

palette4 <- brewer.pal(8,"Dark2")
palette4
# We can use a quick pie chart to see the colors:
pie(rep(1,length(palette1)),col=palette1)
pie(rep(1,length(palette2)),col=palette2)
# We can just stick a few color palettes together
big_palette <- c(palette1, palette2, palette4)
pie(rep(1,length(big_palette)),col=big_palette)

distinctive_palette <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
pie(rep(1,length(distinctive_palette)),col=distinctive_palette)
#Eliminate a couple of similar colours
distinctive_palette <- distinctive_palette[-(21)]
distinctive_palette <- distinctive_palette[-(18)]
distinctive_palette <- distinctive_palette[-(10)]

coul <- c(palette1, palette2)
coul = sample(coul)
pie(rep(1,length(coul)),col=coul)
# Pie chart of all the colors together:
coul = colorRampPalette(c(palette1, palette2, palette3, palette4))(32)
coul2 = sample(coul)
pie(rep(1,length(coul2)),col=coul2)
saved_combo <- coul2

library(pals)
my_glasbey <- glasbey(31)
pie(rep(1,length(my_glasbey)),col=my_glasbey)
#Remember to set GWAS catalog to 1 on each unique figure.
#Retrieve counts of non-unique evidence source for all the SNPs.
all_SNP_non_unique <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% count(study_type)
write.table(all_SNP_non_unique, "all_SNP_non_unique_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_SNP_unique <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% select(study_id, study_type)  %>% distinct() %>% count(study_type)
write.table(all_SNP_unique, "all_SNP_unique_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_gene_non_unique <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% count(study_type)
write.table(all_gene_non_unique, "all_gene_non_unique_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_gene_unique <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% select(study_id, study_type)  %>% distinct() %>% count(study_type)
write.table(all_gene_unique, "all_gene_unique_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

#Do the same for independent loci.
all_SNP_non_unique_per_locus <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>%  group_by(index_SNP_rsid, study_type) %>% summarise(count = n())
write.table(all_SNP_non_unique_per_locus, "all_SNP_non_unique_per_locus_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_SNP_unique_per_locus <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% select(index_SNP_rsid, study_id, study_type) %>% distinct() %>% group_by(index_SNP_rsid, study_type) %>% summarise(length(unique(study_id)))
write.table(all_SNP_unique_per_locus, "all_SNP_unique_per_locus_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_gene_non_unique_per_locus <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>%  group_by(index_SNP_rsid, study_type) %>% summarise(count = n())
write.table(all_gene_non_unique_per_locus, "all_gene_non_unique_per_locus_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

all_gene_unique_per_locus <- all_normal_together_tbl %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% select(index_SNP_rsid, study_id, study_type) %>% distinct() %>% group_by(index_SNP_rsid, study_type) %>% summarise(length(unique(study_id)))
write.table(all_gene_unique_per_locus, "all_gene_unique_per_locus_study_type.txt", quote=F, col.names=T, sep="\t", row.names=F)

#Read in the files above.
all_SNP_non_unique <- read.delim("all_SNP_non_unique_study_type.txt", header=T)
colnames(all_SNP_non_unique)[2] <- "count"
all_SNP_non_unique[all_SNP_non_unique$study_type=='GWAS Catalog',]$count <- 1

all_SNP_unique <- read.delim("all_SNP_unique_study_type.txt", header=T)
colnames(all_SNP_unique)[2] <- "count"
all_SNP_unique[all_SNP_unique$study_type=='GWAS Catalog',]$count <- 1

all_gene_non_unique <- read.delim("all_gene_non_unique_study_type.txt", header=T)
colnames(all_gene_non_unique)[2] <- "count"
all_gene_unique <- read.delim("all_gene_unique_study_type.txt", header=T)
colnames(all_gene_unique)[2] <- "count"

all_SNP_non_unique_per_locus <- read.delim("all_SNP_non_unique_per_locus_study_type.txt", header=T)
all_SNP_non_unique_per_locus[all_SNP_non_unique_per_locus$study_type=='GWAS Catalog',]$count <- 1
all_SNP_unique_per_locus <- read.delim("all_SNP_unique_per_locus_study_type.txt", header=T)
colnames(all_SNP_unique_per_locus)[3] <- "count"
all_SNP_unique_per_locus[all_SNP_unique_per_locus$study_type=='GWAS Catalog',]$count <- 1
all_gene_non_unique_per_locus <- read.delim("all_gene_non_unique_per_locus_study_type.txt", header=T)
all_gene_unique_per_locus <-  read.delim("all_gene_unique_per_locus_study_type.txt", header=T)
colnames(all_gene_unique_per_locus)[3] <- "count"

#Create a fake plot which has all the possible factors in it.
all <- rbind(all_gene_unique, all_SNP_unique)
all_together <- aggregate(count ~ study_type, data=all, sum)
#Create Custom colour palette.
cols <- my_glasbey
names(cols) <- all_together$study_type


#Convert to long format for plotting
plot_all <- function(x, y) {
x$Model <- y
ggplot(x, aes(x=Model, y=count)) + 
  geom_bar(stat="identity", aes(fill = study_type)) + scale_fill_manual(values = cols, drop=TRUE) +  scale_y_continuous(labels = scales::comma) + theme(axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12), legend.title=element_text(face="bold",size=12))
my_pdf <- paste(y, ".pdf", sep="")
my_png <- paste(y, ".png", sep="")
ggsave(my_pdf, plot=last_plot())
ggsave(my_png, plot=last_plot())
}

plot_all(all_SNP_non_unique, 'all_SNP_non_unique')
plot_all(all_SNP_unique, 'all_SNP_unique')
plot_all(all_gene_non_unique, 'all_gene_non_unique')
plot_all(all_gene_unique, 'all_gene_unique')

#Plot each individual locus separately
all_my_index_loci <- unique(all_SNP_non_unique_per_locus$index_SNP_rsid)
for (i in all_my_index_loci)
{
my_out <- paste(i, '_all_SNP_non_unique_per_locus', sep="")
plot_all(all_SNP_non_unique_per_locus[all_SNP_non_unique_per_locus$index_SNP_rsid==i,], my_out)
my_out <- paste(i, '_all_SNP_unique_per_locus', sep="")
plot_all(all_SNP_unique_per_locus, my_out)
my_out <- paste(i, '_all_gene_non_unique_per_locus', sep="")
plot_all(all_gene_non_unique_per_locus, my_out)
my_out <- paste(i, '_all_gene_unique_per_locus', sep="")
plot_all(all_gene_unique_per_locus, my_out)
}


#Do the same for top results.
#Read in the evidence for the top 3 genes.
#Read in the evidence for the top 1 gene.
#Read in top ranked genes for each locus - Model13.
model13_gene <- read.delim("Model13_gene_ranked.txt", sep="\t", header=T)
model13_gene_top1 <- model13_gene[model13_gene$rank==1,]
model13_gene_top3 <- model13_gene[model13_gene$rank < 4,]
#Read in top ranked SNP for each locus - Model09.
model09_snp <- read.delim("Model09_snp_ranked.txt", sep="\t", header=T)
model09_snp_top1 <- model09_snp[model09_snp$rank==1,]
model09_snp_top3 <- model09_snp[model09_snp$rank < 4,]

#Read in top ranked SNP for each locus - Model10.
model10_snp <- read.delim("Model10_snp_ranked.txt", sep="\t", header=T)
model10_snp_top1 <- model10_snp[model10_snp$rank==1,]
model10_snp_top3 <- model10_snp[model10_snp$rank < 4,]

#Read in top ranked SNP for each locus - Model14.
model14_snp <- read.delim("Model14_snp_ranked.txt", sep="\t", header=T)
model14_snp_top1 <- model14_snp[model14_snp$rank==1,]
model14_snp_top3 <- model14_snp[model14_snp$rank < 4,]


loadFile <- function(x) {
  print (x)
  df <- read.delim(x, header=T, stringsAsFactors=F,row.names=NULL, sep="\t")
  df <- tbl_df(df)
}

#Read in the relevant results
model09_snp_top3$file <- paste("../new_study_type/Model09_top_10_snp/", model09_snp_top3$index_SNP_rsid, '_', model09_snp_top3$current_SNP_rsid, ".snps", sep="")
model09_snp_top1$file <- paste("../new_study_type/Model09_top_10_snp/", model09_snp_top1$index_SNP_rsid, '_', model09_snp_top1$current_SNP_rsid, ".snps", sep="")
model09_snp_top3_df <- lapply(model09_snp_top3$file,loadFile)
model09_snp_top1_df <- lapply(model09_snp_top1$file,loadFile)

model10_snp_top3$file <- paste("../new_study_type/Model10_top_10_snp/", model10_snp_top3$index_SNP_rsid, '_', model10_snp_top3$current_SNP_rsid, ".snps", sep="")
model10_snp_top1$file <- paste("../new_study_type/Model10_top_10_snp/", model10_snp_top1$index_SNP_rsid, '_', model10_snp_top1$current_SNP_rsid, ".snps", sep="")
model10_snp_top3_df <- lapply(model10_snp_top3$file,loadFile)
model10_snp_top1_df <- lapply(model10_snp_top1$file,loadFile)

model14_snp_top3$file <- paste("../new_study_type/Model14_top_10_snp/", model14_snp_top3$index_SNP_rsid, '_', model14_snp_top3$current_SNP_rsid, ".snps", sep="")
model14_snp_top1$file <- paste("../new_study_type/Model14_top_10_snp/", model14_snp_top1$index_SNP_rsid, '_', model14_snp_top1$current_SNP_rsid, ".snps", sep="")
model14_snp_top3_df <- lapply(model14_snp_top3$file,loadFile)
model14_snp_top1_df <- lapply(model14_snp_top1$file,loadFile)

model13_gene_top3$file <- paste("../new_study_type/Model13_top_10_gene/", model13_gene_top3$index_SNP_rsid, '_', model13_gene_top3$gene_name, ".genes", sep="")
model13_gene_top1$file <- paste("../new_study_type/Model13_top_10_gene/", model13_gene_top1$index_SNP_rsid, '_', model13_gene_top1$gene_name, ".genes", sep="")
model13_gene_top3_df <- lapply(model13_gene_top3$file,loadFile)
model13_gene_top1_df <- lapply(model13_gene_top1$file,loadFile)

#Process the results and print.
process_print_top <- function(x, y) {
  all_SNP_non_unique <- x %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% dplyr::count(study_type)
  all_SNP_unique <- x %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(current_SNP_rsid)) %>% select(study_id, study_type)  %>% distinct() %>% dplyr::count(study_type) 
  my_locus <- unique(x$index_SNP_rsid)
  my_out <- paste(my_locus, "_", y,  "_unique", sep="")
  colnames(all_SNP_unique)[2] <- 'count'
  plot_all(all_SNP_unique, my_out)
  my_out2 <- paste(my_locus, "_", y, "_non-unique", sep="")
  colnames(all_SNP_non_unique)[2] <- 'count'
  plot_all(all_SNP_non_unique, my_out2)
}

for (a in seq_along(model09_snp_top3_df)) {
  process_print_top(model09_snp_top3_df[[a]], 'model09_snp_top3')  
}

for (a in seq_along(model09_snp_top1_df)) {
  process_print_top(model09_snp_top1_df[[a]], 'model09_snp_top1')  
}

for (a in seq_along(model10_snp_top3_df)) {
  process_print_top(model10_snp_top3_df[[a]], 'model10_snp_top3')   
}
for (a in seq_along(model10_snp_top1_df)) {
  process_print_top(model10_snp_top1_df[[a]], 'model10_snp_top1')   
}

for (a in seq_along(model14_snp_top3_df)) {
  process_print_top(model14_snp_top3_df[[a]], 'model14_snp_top3')   
}
for (a in seq_along(model14_snp_top1_df)) {
  process_print_top(model14_snp_top1_df[[a]], 'model14_snp_top1')   
}


#Process the results and print for genes.
process_print_top2 <- function(x, y) {
  all_SNP_non_unique <- x %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% count(study_type)
  all_SNP_unique <- x %>% filter(!is.na(index_SNP_rsid)) %>% filter(!is.na(gene_name)) %>% select(study_id, study_type)  %>% distinct() %>% count(study_type) 
  my_locus <- unique(x$index_SNP_rsid)
  my_out <- paste(my_locus, "_", y,  "_unique", sep="")
  colnames(all_SNP_unique)[2] <- 'count'
  plot_all(all_SNP_unique, my_out)
  my_out2 <- paste(my_locus, "_", y, "_non-unique", sep="")
  colnames(all_SNP_non_unique)[2] <- 'count'
  plot_all(all_SNP_non_unique, my_out2)
}

for (a in seq_along(model13_gene_top3_df)) {
  process_print_top2(model13_gene_top3_df[[a]], 'model13_gene_top3')  
}

for (a in seq_along(model13_gene_top1_df)) {
  process_print_top2(model13_gene_top1_df[[a]], 'model13_gene_top1')  
}
