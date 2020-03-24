library("tidyverse")
library("extrafont")
font_import()
loadfonts()
remotes::install_github('noamross/viridisLite@turbo')

#List of loci included in the paper
paper_loci <- c("1q21.3a", "1q21.3b", "2p13.3", "2q12.1", "2q37.1", "4q27", "5p13.2", "5q31.1", "6p21.32", "6p21.33", "8q21.13", "10p15.1", "10q21.2", "11p13", "11q13.1", "11q13.5", "11q24.3", "12q15", "14q13.2", "14q32.32", "16p13.13", "17q21.2", "17q25.3", "19p13.2", "20q13.33")
#read in file Gene ranking
gene_ranking <- read.delim("Model13_gene_ranked.txt", header=T, stringsAsFactors = FALSE)
#read in file with locus-top SNP mapping
locus_names <- read.delim("Paternoster2015_locus_names.txt", header=T, sep="\t", stringsAsFactors = FALSE)
#read in a file with lead SNPs to keep in the table (taken from Table 1)
lead_SNPs <- read.delim("index_SNPs_to_keep.txt", header=T, sep="\t", stringsAsFactors=FALSE)
#read in coloc effect direction
coloc_effect <- read.delim("coloc_effect_loci.txt", header=T, sep="\t", stringsAsFactors = FALSE)
#Merge loci with gene ranking/scores
gene_ranking_merged <- merge(gene_ranking, locus_names, by.x="index_SNP_rsid", by.y="Top.SNP", all.x=T)
#Filter based on rank (>4) and loci to be included in the paper
gene_ranking_merged <- gene_ranking_merged[-c(6)]
gene_ranking_merged <- gene_ranking_merged[gene_ranking_merged$rank < 4,]
gene_ranking_merged <- gene_ranking_merged[gene_ranking_merged$Locus %in% paper_loci,]



loadFile <- function(x) {
  gene_name <- x["gene_name"]
  index_SNP <- x["index_SNP_rsid"]
  locus <- x["Locus"]
  data_file <- paste("old_study_type/Model13_top_10_gene/",index_SNP,"_", gene_name, ".genes", sep="")
  df <- read.delim(data_file, header=T, stringsAsFactors=F, row.names=NULL, sep="\t")
  df
}

all_filtered <- apply(gene_ranking_merged, 1, loadFile)
all_filtered_together <- do.call(rbind,all_filtered)
all_filtered_together <- as.data.frame(all_filtered_together)
#Recode to more general categories used in generating scores
all_filtered_together <- tbl_df(all_filtered_together)
#all_filtered_together <- all_filtered_together %>% mutate(study_type = recode(study_type, 'oe' = 'enhancer', 'promoter enhancer'='promoter-enhancer','y'='enhancer','enhancers'='enhancer','bait'='promoter', 'pir'='enhancer','x'='promoter','regulatory variant predictions'='regulatory variant prediction', 'enhancer-promoter'='promoter-enhancer', 'coloc - SNP'='coloc', 'FAIRE'='FAIRE-seq', 'DEG'='DGE', 'CTCF sites'='CTCF', 'coloc - gene'='coloc', 'network analysis'='PrixFixe', 'ChIP-seq peaks'='ChIP-Seq', 'epigenome'='active chromatin state', 'Converved genomic regulatory blocks'='Conserved genomic regulatory blocks', 'ChIP-Seq'='ChIP-seq'))

all_filtered_together <- all_filtered_together %>% mutate(study_type = recode(study_type, 'oe' = 'promoter-enhancer', 'MTGDR classifier' = 'DGE', 'promoter enhancer'='promoter-enhancer', 'enhancer'='promoter-enhancer','y'='promoter-enhancer','enhancers'='promoter-enhancer','bait'='promoter-enhancer', 'DGE meta'='DGE', 'pir'='promoter-enhancer','x'='promoter-enhancer','regulatory variant predictions'='regulatory variant prediction', 'enhancer-promoter'='promoter-enhancer', 'sQTL'='eQTL', 'coloc - SNP'='coloc', 'exon-ratio eQTL'='eQTL', 'polyA ratio level eQTL'='eQTL', 'FAIRE'='FAIRE-seq', 'isoQTL'='eQTL', 'DEG'='DGE', 'CTCF sites'='CTCF', 'aseQTL'='eQTL', 'coloc - gene'='coloc', 'exon eQTL'='eQTL', 'context specific eQTL'='eQTL', 'ase hQTL'='hQTL', 'transcript ratio eQTL'='eQTL', 'network analysis'='PrixFixe', 'DNMT3B_peaks' = "ChIP-seq", 'ChIP-seq peaks'='ChIP-Seq', 'DNMT3A_peaks'='ChIP-Seq',  'epigenome'='active chromatin state', 'miRNA target sites'='miRNA', 'miRNA seed regions'='miRNA','ChIP-Seq'='ChIP-seq', 'Converved genomic regulatory blocks'='Conserved genomic regulatory blocks', 'ChIP-Seq'='ChIP-seq'))
all_study_types <- unique(all_filtered_together$study_type)

#Sum scores in each category - not needed for the figure below       
sum_final_scores <- function(x) {
#Filter down to entries with a given study_type
  study_type_scores <- all_filtered_together[all_filtered_together$study_type==x,]$final_score
  sum_of_scores <- sum(study_type_scores)
  sum_of_scores
}
summed_scores <- sapply(all_study_types, sum_final_scores)
sort(summed_scores, decreasing=T)

#Adjustment to do with study variety
adjustment <- all_filtered_together %>% group_by(index_SNP_rsid, gene_name) %>% summarise(sqrt((length(unique(study_type)) + length(unique(study_id)))/2))
names(adjustment)[3] <- "adjustment_multiplier"

#Sum TWAS coloc results with TWAS, and within other cats. Rename promoter-enhancer to HiC
total_scores <- all_filtered_together %>% mutate(study_type = recode(study_type, 'TWAS coloc' = 'TWAS', "promoter-enhancer" = 'Hi-C')) %>% group_by(index_SNP_rsid, gene_name, study_type) %>% summarise(sum(final_score)) 
names(total_scores)[4] <- "study_type_score"

#Merge adjustment to the summed scores for each cat to generate final scores for each cat
merge_results <- merge(total_scores, adjustment, by=c("index_SNP_rsid", "gene_name"))
merge_results$final_study_type_score <- merge_results$adjustment_multiplier * merge_results$study_type_score


top_10_categories <- c("coloc", "TWAS", "Hi-C", "eQTL", "DGE", "regfm", "PrixFixe", "mQTL", "pQTL", "hQTL")

#Filter results to include only the top 10 categories
merge_results_top10 <- merge_results[merge_results$study_type %in% top_10_categories,]
#And further filtering based on rsIDs -  select lead SNPs, same as Table 1
merge_results_top10 <- merge_results_top10[merge_results_top10$index_SNP_rsid %in% lead_SNPs$index_SNP_rsid,]

#Drop unneccessary columns
final_df_top10 <- tbl_df(na.omit(merge_results_top10[-c(4, 5)]))

#Convert to wide format
final_df_top10_wide <- spread(final_df_top10, study_type, final_study_type_score)

#Add final_score (total) for each gene
final_df_top10_wide2 <- merge(final_df_top10_wide,gene_ranking_merged, by=c("index_SNP_rsid", "gene_name"))

#Reorder columns
final_df_top10_wide3 <- final_df_top10_wide2[c("rank", "Locus", "index_SNP_rsid", "gene_name", "coloc", "TWAS", "eQTL", "Hi-C", "DGE", "regfm", "PrixFixe", "mQTL", "pQTL", "hQTL", "final_score")]
#Create a fake temp column for sorting:
final_df_top10_wide3$sorting <- paste(final_df_top10_wide3$index_SNP_rsid, final_df_top10_wide3$rank, sep="_")
#Create a vector to sort by
add_rank <- function (x) {
  return_value = c(paste(x, 1, sep="_"), paste(x, 2, sep="_"), paste(x, 3, sep="_"))
}
my_lists <- lapply(lead_SNPs$index_SNP_rsid, add_rank)
sorting_order <-  unlist(my_lists)
 
#Sort df by vector                          
data <- final_df_top10_wide3[match(sorting_order, final_df_top10_wide3$sorting),]
#Rename final score to total score (per gene)
names(data)[names(data) == "final_score"] <- "total score"

#Drop unneccesary columns and join the other ones into a single identifer.
data$identifier <- paste(data$Locus, " / ", data$index_SNP_rsid, " ", data$gene_name)
rownames(data) <- data$identifier

#Remove also final score (index 15)
data_final <- data[(-c(1:4,15, 16,17))]

#Remove the MHC region
#data <- data[data$Locus != "6p21.32", ]

#Transform to long format again for plotting. Change Id and category to factor so that display order is maintained in the plot
long_data <- data_final
long_data$id <- rownames(long_data)
long_data$id <- factor(long_data$id, levels = long_data$id)
long_data <- tbl_df(long_data)
long_data_final <- long_data %>% gather(key = "category", value = "score", -id)
long_data_final$category <- factor(long_data_final$category, levels=top_10_categories)

#max_value on the graph
1186.02464

#round scores down to integer
long_data_final$score <- round(long_data_final$score,0)
ggplot(long_data_final, aes(category, fct_rev(as_factor(id)), col = score, fill = score, label=score)) +
geom_tile() +
geom_text(col = "white") +
theme_minimal() + theme(axis.title=element_blank(), axis.text.x=element_text(size=rel(2), colour="maroon"), legend.title=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1.5)), legend.key.size=unit(1.75, "cm"), axis.ticks=element_blank()) +
scale_fill_viridis_c(option = "turbo", na.value="white",breaks = c(0,100,200,300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200)) +
scale_color_viridis_c(option = "turbo", na.value="white",breaks = c(0,100,200,300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200))

scale_fill_distiller(palette = "Spectral") + 
scale_color_distiller(palette = "Spectral")
                                                                                                                                                                                               