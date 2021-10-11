library("tidyverse")
library("extrafont")
library("cowplot")
library("hrbrthemes")
library(scales)
library("gtools")
font_import()
loadfonts()

my_input <- tbl_df(read.csv("figure_summary_table.csv", header = TRUE, stringsAsFactors = FALSE))
my_input <- my_input[gtools::mixedorder(my_input$cytoband), ]
#Croup the variables into categories
my_bubbles <- c("rsid", "cytoband", "HGNC_symbol", "open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank")
my_heatmap <- c("rsid", "cytoband", "HGNC_symbol", "total_evidence_sources", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome")
my_tick <- c("rsid", "cytoband", "HGNC_symbol", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP")
my_all_ranked <- c("rsid", "cytoband", "HGNC_symbol",  "open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP")
my_input_bubbles <- my_input[,my_bubbles]
my_input_heatmap <- my_input[,my_heatmap]
my_input_tick <- my_input[,my_tick]
my_input_ranked <- my_input[,my_all_ranked]

keycol <- "method" 
valuecol <- "count"
gathercols <- c("total_evidence_sources", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome")
heatmap_long <- gather_(my_input_heatmap, keycol, valuecol, gathercols)
heatmap_long$id <- paste(heatmap_long$cytoband, " / ", heatmap_long$rsid,  " / ", heatmap_long$HGNC_symbol)
heatmap_long$id <- factor(heatmap_long$id, levels = unique(heatmap_long$id))
heatmap_long$method <- factor(heatmap_long$method, levels=unique(heatmap_long$method))
heatmap_long$count_no_0 <- gsub("\\<0\\>", "", heatmap_long$count)

#Dimensions 6 X 7 inches.
#Heatmap with expression resources containing also the total scores.
heatmap_all <- ggplot(heatmap_long, aes(method, fct_rev(as_factor(id)), col = count, fill = count, label=count_no_0)) +
geom_tile(color = "gray") +
geom_text(col = "white", size = 2) +
theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="maroon", angle = 90, 
hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), legend.key.size=unit(15, "pt"),
axis.ticks=element_blank()) +  
scale_x_discrete(breaks=c("total_evidence_sources", "total_evidence_pieces", "coloc", "smultixcan",  "smr", "dge_gxp", "dge_proteome"), 
labels=c("Total evidence sources", "Total evidence pieces", "coloc", "SMultiXcan",  "SMR", "transcriptome", "proteome"), position = "top") +
scale_color_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
scale_fill_viridis_c(option = "turbo", na.value="white", breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
ggsave("heatmap_score_all.pdf", heatmap_all, dpi=300, height=7, width=6, units="in")


#Bubble plot for ranks and binary
#Dimensions 6 X 7 inches
keycol <- "method" 
valuecol <- "rank"
gathercols <- c("open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP")
non_ranked <- c( "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP")
ranked_long <- gather_(my_input_ranked, keycol, valuecol, gathercols)
ranked_long$id <- paste(ranked_long$cytoband, " / ", ranked_long$rsid,  " / ", ranked_long$HGNC_symbol)
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
scale_x_discrete(breaks=c("open_targets_prioritization_rank", "pops_prioritization_rank", "POSTGAP_prioritization_rank", "DEPICT_prioritization", "MAGMA_prioritization", "MendelVar_sig_enrichment", "MendelVar_skin_keywords", "VEP"), 
labels=c("Open Targets", "PoPs", "POSTGAP", "DEPICT", "MAGMA", "MendelVar enrichment", "MendelVar skin", "VEP"), 
position = "top")  +
scale_size(range = c(4, 2), breaks=c(1, 2, 3), name="Rank") +
guides(shape=FALSE, color=FALSE) 
ggsave("bubble_plot.pdf", bubble, dpi=300, height=7, width=6, units="in")


#Barchart for score and number of evidence
#Output 8 x 7 inches
barchart_in <- heatmap_long[heatmap_long$method %in% c("total_evidence_sources", "total_evidence_pieces"), ]
barchart_in$method <-gsub("total_evidence_sources", "Total evidence sources", barchart_in$method )
barchart_in$method <-gsub("total_evidence_pieces", "Total evidence pieces", barchart_in$method )
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
ggsave("barchart_score.pdf", barchart_all, dpi=300, height=7, width=8, units="in")

##Heatmap without evidence score - expression-based resources
heatmap_short <- heatmap_long %>% filter(!(method %in% c("total_evidence_sources", "total_evidence_pieces")))
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
ggsave("heatmap_score_select.pdf", heatmap_select, dpi=300, height=7, width=6, units="in")


