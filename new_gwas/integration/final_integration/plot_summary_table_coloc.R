library(tools) 
library("tidyverse")
library(viridis)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least 1 argument must be supplied", call.=FALSE)}

coloc_results <- args[1]
cyto_file <- args[2] 
gwas_name <- args[3]

coloc_t <- 0.95

my_input <- read.delim(coloc_results, header = TRUE, stringsAsFactors = FALSE)
#Provide the original file with rsid to cytoband mapping
my_cyto <- read.delim(cyto_file, stringsAsFactors = F, header=F)
my_cyto$cytoband <- paste(my_cyto$V1, my_cyto$V8, sep="")
my_cyto = subset(my_cyto, select = -c(V1, V2, V3, V5, V6, V7, V8, V9))
colnames(my_cyto)[1] <- "rsid"
my_input <- merge(my_input, my_cyto, by="rsid")

#Filter to keep only results above a certain H4 threshold#
select <- my_input[my_input$PP.H4.abf > coloc_t,]
select <- select[gtools::mixedorder(select$cytoband), ]
select$id <- paste(select$cytoband, " / ", select$rsid,  " / ", select$hugo_name)
#Create labels for conditions
select$label <- paste(select$tissue, select$condition, sep=" ")
select$label <- gsub(" NA", "", select$label, ignore.case = FALSE)
select$label <- gsub("_", " ", select$label, ignore.case = FALSE)
select$label <- gsub(" naive", "", select$label, ignore.case = FALSE)
select$label <- gsub("IFNg IFNg", "IFNg", select$label, ignore.case = FALSE)
select$label <- gsub("B-Cell CD19", "B-Cell", select$label, ignore.case = FALSE)
select$label <- toTitleCase(select$label)
select$label <- gsub("Anti-CD3-CD28 Anti-CD3-CD28", "Anti-CD3-CD28", select$label, ignore.case = FALSE)
select$label <- gsub("Cells EBV-Transformed Lymphocytes", "LCL", select$label, ignore.case = FALSE)
select$label <- gsub("Cells Cultured Fibroblasts", "Fibroblast", select$label, ignore.case = FALSE)
select$label <- gsub("IFNg\\+Salmonella", "", select$label, ignore.case = FALSE)
select$label <- gsub("Salmonella Salmonella", "Salmonella", select$label, ignore.case = FALSE)
select$label <- gsub("Listeria Listeria", "Listeria", select$label, ignore.case = FALSE)
select$label <- gsub("IFN24 IFNg", "IFNg", select$label, ignore.case = FALSE)
select$label <- gsub("LPS2 LPS", "LPS", select$label, ignore.case = FALSE)
select$label <- gsub("R848 R848", "R848", select$label, ignore.case = FALSE)
select$label <- gsub("LPS LPS", "LPS", select$label, ignore.case = FALSE)
select$label <- gsub("Pam3CSK4 Pam3CSK4", "Pam3CSK4", select$label, ignore.case = FALSE)
select$label <- gsub("CD4 T-Cell", "T-Cell CD4", select$label, ignore.case = FALSE)
select$label <- gsub("CD8 T-Cell", "T-Cell CD8", select$label, ignore.case = FALSE)
select$label <- gsub("Transverse Colon", "Colon Transverse", select$label, ignore.case = FALSE)
select$label <- gsub("Small Intestine Terminal Ileum", "Ileum", select$label, ignore.case = FALSE)
select$label <- gsub("Skin_Sun_Exposed_Lower_leg", "Sun-exposed skin", select$label, ignore.case = FALSE)
select$label <- gsub("Skin not Sun Exposed Suprapubic", "Unexposed skin", select$label, ignore.case = FALSE)
select$study <- gsub("_", "", select$study)
select$detailed <- paste(select$label, " (", select$study, ")", sep="")
select <- as_tibble(select)
#Subset to different tissue types.
#Whole blood 
select_blood <- select %>% filter(label %in% c("Whole Blood", "Blood", "plasma"))

#Immune cell types
select_immune <- select %>% filter(study %in% c("Alasoo_2018", "BLUEPRINT", "Fairfax_2012", "Fairfax_2014", "GEUVADIS", "Kasela_2017", "Naranbhai_2015", "Nedelec_2016", "Quach_2016", "Schmiedel_2018") | tissue %in% c("LCL", "T-cell_CD8", "T-cell_CD4", "monocyte_CD14", "neutrophil_CD15", "B-cell_CD19", "Cells_EBV-transformed_lymphocytes", "T-cell"))

#Skin and related
select_skin <- select %>% filter(tissue %in% c("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic", "Cells_Cultured_fibroblasts", "skin", "Fibroblast"))

#Immune-related organs and tissues
select_organs <- select %>% filter(tissue %in% c("Lung", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Mucosa", "Small_Intestine_Terminal_Ileum", "Spleen", "transverse_colon", "rectum", "ileum"))

#Create labels for each category and order alphabetically
select_blood$label <- "Whole Blood"
my_skin_label <- sort(unique(select_skin$label))
my_immune_label <- sort(unique(select_immune$label))
my_organ_label <- sort(unique(select_organs$label))

#Bubble plot
#Blood
blood_all <- ggplot(select_blood, aes(x=study, y=fct_rev(as_factor(id)), color=PP.H4.abf)) +
geom_point(size=6) +
#theme_minimal() + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="red", 
angle = 90, hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), 
legend.key.size=unit(15, "pt"), 
#axis.ticks.length=unit(0.5,"cm"),
#axis.ticks=element_blank(),
#axis.line.x = element_line(color="black", size = 0.2),
#axis.line.y = element_line(color="black", size = 0.2),
panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
labs(color="PPH4")+ 
scale_color_viridis(option = "D", limits=c(coloc_t, 1))
figure_output0 <- paste("coloc_blood_all_", gwas_name, ".pdf", sep="")
ggsave(figure_output0, blood_all, dpi=300, height=25, width=6, units="in")

#Skin
skin_all <- ggplot(select_skin, aes(x=detailed, y=fct_rev(as_factor(id)), color=PP.H4.abf)) +
geom_point(size=6) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="red", 
angle = 90, hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), 
legend.key.size=unit(15, "pt"), 
panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
labs(color="PPH4")+ 
scale_color_viridis(option = "D", limits=c(coloc_t, 1))

figure_output <- paste("coloc_skin_all_", gwas_name, ".pdf", sep="")
ggsave(figure_output, skin_all, dpi=300, height=30, width=6, units="in")

##Immune cell types
immune_all <- ggplot(select_immune, aes(x=detailed, y=fct_rev(as_factor(id)), color=PP.H4.abf)) +
geom_point(size=6) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="red", 
angle = 90, hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), 
legend.key.size=unit(15, "pt"), 
panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
labs(color="PPH4")+ 
scale_color_viridis(option = "D", limits=c(coloc_t, 1))
figure_output2 <- paste("coloc_immune_all_", gwas_name, ".pdf", sep="")
ggsave(figure_output2, immune_all, dpi=300, height=40, width=6, units="in")

#Organs
organs_all <- ggplot(select_organs, aes(x=detailed, y=fct_rev(as_factor(id)), color=PP.H4.abf)) +
geom_point(size=6) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.title=element_blank(), axis.text.x=element_text(colour="red", 
angle = 90, hjust = 0), legend.title=element_text(size = 10), legend.text=element_text(size=8), 
legend.key.size=unit(15, "pt"), 
panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
labs(color="PPH4")+ 
scale_color_viridis(option = "D", limits=c(coloc_t, 1))
figure_output3 <- paste("coloc_organs_all_", gwas_name, ".pdf", sep="")
ggsave(figure_output3, organs_all, dpi=300, height=35, width=6, units="in")
