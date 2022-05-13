library("tidyverse")
library("extrafont")
library("cowplot")
#library("hrbrthemes")
library(scales)
library("gtools")
library("tools")
font_import()
loadfonts()

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE)}

gwas_name <- args[1]
#File classifying SNPs as known and novel
known_novel <- args[2]
known_novel_df <- read.delim(known_novel, stringsAsFactors=F)
known_loci <- known_novel_df[known_novel_df$Known...Novel=="Known",]
novel_loci <- known_novel_df[known_novel_df$Known...Novel=="Novel",]

input_file <- paste("figure_summary_table_", gwas_name, ".csv", sep="")
my_input <- tbl_df(read.csv(input_file, header = TRUE, stringsAsFactors = FALSE))
my_input <- my_input[gtools::mixedorder(my_input$cytoband), ]
#Croup the variables into categories
my_bubbles <- c("rsid", "cytoband", "HGNC_symbol", "CHR", "POS","open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank")
my_heatmap <- c("rsid", "cytoband", "HGNC_symbol", "CHR", "POS", "total_evidence_types", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome")
my_tick <- c("rsid", "cytoband", "HGNC_symbol", "CHR", "POS", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP_intron", "VEP_missense")
my_all_ranked <- c("rsid", "cytoband", "HGNC_symbol", "CHR", "POS", "open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP_intron", "VEP_missense")
my_input_bubbles <- my_input[,my_bubbles]
my_input_heatmap <- my_input[,my_heatmap]
my_input_tick <- my_input[,my_tick]
my_input_ranked <- my_input[,my_all_ranked]

keycol <- "method" 
valuecol <- "count"
gathercols <- c("total_evidence_types", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome")
heatmap_long <- gather_(my_input_heatmap, keycol, valuecol, gathercols)
heatmap_long$id <- paste(heatmap_long$cytoband, " / ", heatmap_long$rsid,  " / ", heatmap_long$HGNC_symbol)
#Sort by chromosomal position
heatmap_long <- arrange(heatmap_long, CHR, POS)
heatmap_long$id <- factor(heatmap_long$id, levels = unique(heatmap_long$id))
heatmap_long$method <- factor(heatmap_long$method, levels=unique(heatmap_long$method))
heatmap_long$count_no_0 <- gsub("\\<0\\>", "", heatmap_long$count)

#Dimensions 6 X 7 inches.
#Heatmap with expression retypes containing also the total scores.
heatmap_all <- ggplot(heatmap_long, aes(method, fct_rev(as_factor(id)), col = count, fill = count, label=count_no_0)) +
  geom_tile(color = "gray") +
  geom_text(col = "white", size = 2) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="maroon", angle = 90, 
                                                                                                                   hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), legend.key.size=unit(15, "pt"),
                          axis.ticks=element_blank()) +  
  scale_x_discrete(breaks=c("total_evidence_types", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome"), 
                   labels=c("Evidence types", "Total evidence score", "coloc", "SMultiXcan",  "SMR", "transcriptome", "proteome"), position = "top") +
  scale_color_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_fill_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
figure_output <- paste("heatmap_score_all_", gwas_name, ".pdf", sep="")
ggsave(figure_output, heatmap_all, dpi=300, height=15, width=6, units="in")


#Bubble plot for ranks and binary
#Dimensions 6 X 7 inches
keycol <- "method" 
valuecol <- "rank"
gathercols <- c("open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP_intron", "VEP_missense")
non_ranked <- c( "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP_intron", "VEP_missense")
ranked_long <- gather_(my_input_ranked, keycol, valuecol, gathercols)
ranked_long$id <- paste(ranked_long$cytoband, " / ", ranked_long$rsid,  " / ", ranked_long$HGNC_symbol)
#Sort by chromosomal position
ranked_long <- arrange(ranked_long, CHR, POS)
ranked_long$id <- factor(ranked_long$id, levels = unique(ranked_long$id))
ranked_long$method <- factor(ranked_long$method, levels=unique(ranked_long$method))
ranked_long$rank <- ifelse(ranked_long$rank != 0, ranked_long$rank , NA)
ranked_long[ranked_long$method %in% non_ranked,]$rank <- ifelse(ranked_long[ranked_long$method %in% non_ranked,]$rank == 1, 2 , ranked_long[ranked_long$method %in% non_ranked,]$rank)
ranked_long$group <- ifelse(ranked_long$method %in% non_ranked, "a" , "b")

bubble <- ggplot(ranked_long, aes(method, fct_rev(as_factor(id)), size = rank, label=rank, shape=group, color=group)) +
  geom_point() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="maroon", 
                                                                                                 angle = 90, hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), 
        legend.key.size=unit(15, "pt"), axis.ticks=element_blank()) + 
  scale_shape_manual(values=c(18, 16)) +
  scale_color_manual(values=c("lightskyblue", "black")) +
  scale_x_discrete(breaks=c("open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP_intron", "VEP_missense"), 
                   labels=c("Open Targets", "PoPs", "POSTGAP", "DEPICT", "MAGMA", "MendelVar enrichment", "MendelVar skin", "VEP intron", "VEP missense"), 
                   position = "top")  +
  scale_size(range = c(4, 2), breaks=c(1, 2, 3), name="Rank") +
  guides(shape=FALSE, color=FALSE) 
figure_output2 <- paste("bubble_plot_", gwas_name, ".pdf", sep="")
ggsave(figure_output2, bubble, dpi=300, height=15, width=6, units="in")

#Barchart for score and number of evidence
#Output 8 x 7 inches
barchart_in <- heatmap_long[heatmap_long$method %in% c("total_evidence_types", "total_evidence_pieces"), ]
barchart_in$method <-gsub("total_evidence_types", "Evidence types", barchart_in$method )
barchart_in$method <-gsub("total_evidence_pieces", "Total evidence score", barchart_in$method )
#barchart_in <- barchart_in %>% map_df(rev)

barchart_all <- ggplot(barchart_in, aes(x = reorder(id, desc(id)), y = count)) +
  geom_col(aes(fill = method)) +
  facet_wrap(~ method) + 
  coord_flip() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title=element_blank(), legend.position = "none", 
        axis.ticks.y = element_blank(), axis.line.x = element_line(color="black", size = 0.2)) + 
  scale_y_continuous( breaks=pretty_breaks(n=6)) + 
  scale_x_discrete(limits = rev(levels(id)))
figure_output3 <- paste("barchart_score_", gwas_name, ".pdf", sep="")
ggsave(figure_output3, barchart_all, dpi=300, height=15, width=8, units="in")

#Subset the barchart to known and novel loci
barchart_known <- barchart_in[barchart_in$rsid %in% known_loci$RSID,]
barchart_novel <- barchart_in[barchart_in$rsid %in% novel_loci$RSID,]
#Plot Known and novel loci seperately. Plot total evidence score and type seperately too.
barchart_known_score <- barchart_known[barchart_known$method == "Total evidence score",]
barchart_known_type <- barchart_known[barchart_known$method == "Evidence types",]
barchart_novel_score <- barchart_novel[barchart_novel$method == "Total evidence score",]
barchart_novel_type <- barchart_novel[barchart_novel$method == "Evidence types",]

barchart_known_score_fig <- ggplot(barchart_known_score, aes(x = reorder(id, desc(id)), y = count)) +
  geom_col(aes(fill = "peru")) + 
  coord_flip() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title=element_blank(), legend.position = "none", 
        axis.ticks.y = element_blank(), axis.line.x = element_line(color="black", size = 0.2)) + 
  scale_y_continuous( breaks=pretty_breaks(n=6)) + 
  scale_x_discrete(limits = rev(levels(id)))
figure_output3 <- paste("barchart_score_known_", gwas_name, ".pdf", sep="")
ggsave(figure_output3, barchart_known_score_fig, dpi=300, height=12, width=8, units="in")

##Heatmap without evidence score - expression-based retypes
heatmap_short <- heatmap_long %>% filter(!(method %in% c("total_evidence_types", "total_evidence_pieces")))
heatmap_select <- ggplot(heatmap_short, aes(method, fct_rev(as_factor(id)), col = count, fill = count, label=count_no_0)) +
  geom_tile(color = "gray") +
  geom_text(col = "white", size = 2) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="maroon", angle = 90, 
                                                                                                                   hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), legend.key.size=unit(15, "pt"),
                          axis.ticks=element_blank()) +  
  scale_x_discrete(breaks=c("coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome"), 
                   labels=c("coloc", "SMultiXcan",  "SMR", "transcriptome", "proteome"), position = "top") +
  scale_color_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_fill_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
figure_output5 <- paste("heatmap_score_select_", gwas_name, ".pdf", sep="")
ggsave(figure_output5, heatmap_select, dpi=300, height=15, width=5, units="in")


