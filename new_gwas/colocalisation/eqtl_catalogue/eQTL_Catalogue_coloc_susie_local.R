library("dplyr")
library("readr")
library("coloc")
#library("GenomicRanges")
#library("Rsamtools")
library("hablar")
library("tools")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least one argument must be supplied (rsid).n", call.=FALSE)}


my_rsid <- args[1]


#This script contains all coloc runs for a given SNP and accepts SNP rsid as an argument.
#The main result (H4 result) also saved in aggregate column - one table per RSID.
#Use all genes within a 1Mbp of each SNP to match settings used in eQTL detection in eQTL Catalog.
#USe imported GTEx as version 8

#Load required functions
source("/newhome/qh18484/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue/eQTL_Catalogue_helper_functions_susie.R")

#Create an empty data frame to store all the colocalisation results.
coloc_all <- data.frame(nsnps=double(),
                        PP.H0.abf = double(), 
                        PP.H1.abf = double(),
                        PP.H2.abf = double(), 
                        PP.H3.abf = double(), 
                        PP.H4.abf = double(), 
                        hugo_name = character(),
                        ensembl_id = character(),
                        rsid = character(),
                        study = character(),
                        tissue = character(),
                        condition = character(),
                        feature = character(),
                        stringsAsFactors=FALSE) 

missing_runs <- data.frame(rsid= character(),
                           ensembl_id = character(),
                           hugo_name = character(),
                           study = character(),
                           tissue_label = character(),
                           condition_label = character(),
                           quant_method = character(),
                           stringsAsFactors=FALSE
  
)

coloc_all <-  dplyr::as_tibble(coloc_all)

#Read in the file with target SNPs and interval ranges lifted over to GRCh38.
my_ranges <- read.table("interval_r2_0.2_1k_hg38.bed", stringsAsFactors = F, header=F)
colnames(my_ranges) <- c("chrom", "start", "end", "rsid")

#Read in the file with 1Mbp interval around index SNPs.
#my_ranges <- read.table("paternoster_2015_index_snps_sorted_1Mbp.bed", stringsAsFactors = F, header=F)
#colnames(my_ranges) <- c("chrom", "start", "end", "rsid")

#Read in file with the list of genes within 1 Mbp of each index SNP.
my_genes <- read.table("paternoster_2015_index_snps_sorted_1Mbp_genes_processed.bed", stringsAsFactors = F, header=F)
colnames(my_genes) <- c("chrom", "start", "end", "rsid", "gene_chrom", "gene_start", "gene_end", "strand", "type", "Ensembl_transcript_id", "Ensembl_gene_id", "hugo_gene_name", "huge_gene_id", "file")

#Lookup my LD interval range. Need to use a file with GRCh38 coordinates.
my_chrom = my_ranges[my_ranges$rsid == my_rsid,]$chrom
my_start = my_ranges[my_ranges$rsid == my_rsid,]$start
my_end= my_ranges[my_ranges$rsid == my_rsid,]$end
#region <- GRanges(my_chrom, IRanges(start=my_start, end=my_end))
coordinates <- paste(my_chrom, ":", my_start, "-", my_end, sep="")

#Lookup my Ensembl gene IDs
my_ensembl_genes = unique(my_genes[my_genes$rsid == my_rsid,]$Ensembl_gene_id)

#Load file with my eczema summary stats for a given variant.
eczema_file = paste("GWAS_intervals/", my_rsid, "_r2_0.2_1k.gwas",sep="")
my_eczema <- read.table(eczema_file, stringsAsFactors = F, header=T)
my_eczema <- my_eczema %>% dplyr::as_tibble() %>% dplyr::rename(MAF = maf) %>%
  hablar::convert(num(CHR, POS, N, MAF, BETA, SE, Z_SCORE, PVAL))

#Load list of eQTL Catalog resources
tabix_paths_0 = read.delim("tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
imported_tabix_paths_0 = read.delim("tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

#Filter both tables to keep only studies and tissues we are interested in.
tabix_paths <- tabix_paths_0 %>% filter(study == "Alasoo_2018" | study == "BLUEPRINT" | study == "CEDAR"| 
study == "Fairfax_2012" | study == "Fairfax_2014" | study == "GENCORD" | study == "GEUVADIS" |
study == "Kasela_2017" | study == "Lepik_2017" | study == "Naranbhai_2015" | study == "Nedelec_2016" |
study == "Quach_2016" | study == "Schmiedel_2018" | (study == "TwinsUK" & qtl_group != "fat"))
imported_tabix_paths <- imported_tabix_paths_0 %>% filter(study == "GTEx_V8" & (qtl_group == "Cells_Cultured_fibroblasts" |
qtl_group == "Cells_EBV-transformed_lymphocytes" | qtl_group == "Colon_Sigmoid" | qtl_group == "Colon_Transverse" |
qtl_group == "Esophagus_Mucosa" | qtl_group == "Lung" | qtl_group ==  "Skin_Not_Sun_Exposed_Suprapubic" |
qtl_group == "Skin_Sun_Exposed_Lower_leg" | qtl_group == "Small_Intestine_Terminal_Ileum" |
qtl_group == "Spleen" | qtl_group == "Whole_Blood"))

#Merge the two tables
tabix_paths_all <- bind_rows(tabix_paths, imported_tabix_paths)

#Not interested in txtrev
tabix_paths_all <- tabix_paths_all %>% filter (quant_method == "ge" | quant_method == "microarray")

#Change file paths to local
tabix_paths_all$ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/\\w+/ge", "catalog", tabix_paths_all$ftp_path)
tabix_paths_all$ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/\\w+/microarray", "catalog", tabix_paths_all$ftp_path)

#Extract column names from first file

#column_names = colnames(readr::read_tsv("catalog/Alasoo_2018_ge_macrophage_naive.all.tsv.gz", n_max = 1))
column_names = colnames(readr::read_tsv("catalog/header.tsv", n_max = 1))
column_names_gtex = colnames(readr::read_tsv("catalog/header_gtex.tsv", n_max = 1))


for (a in 1:nrow(tabix_paths_all)) {
 row <- tabix_paths_all[a,] 
 for (my_gene in my_ensembl_genes) {

#summary_stats <- import_eQTLCatalogue(row$ftp_path, region, selected_gene_id = my_gene, column_names)
if (row$study == "GTEx_V8") {
  summary_stats <- import_eQTLCatalogue_ver_tabix(row$ftp_path, coordinates, selected_gene_id = my_gene, column_names_gtex)
  } 
else {
    summary_stats <- import_eQTLCatalogue_ver_tabix(row$ftp_path, coordinates, selected_gene_id = my_gene, column_names)
}

#Check if df null, If so add to a file containing missing runs. 
my_hugo_name = my_genes[my_genes$Ensembl_gene_id == my_gene,]$hugo_gene_name[1]
if (is.null(summary_stats))
  {
  missing_runs[nrow(missing_runs) + 1,] = c(my_rsid, my_gene, my_hugo_name, row$study, row$tissue_label, row$condition_label, row$quant_method)
  next
  }

#Check if no rows for gene found
if (dim(summary_stats)[1] == 0)
  {
  print ("Feature not present in the eQTL catalogue.")
  next
  }

#Summary stats without rsids
#summary_stats = import_eQTLCatalogue_no_rsid(platelet_df$ftp_path, region, selected_gene_id = my_gene, column_names)

#Merge the GWAS and eQTL file by chrom and pos
#merged <- merge(summary_stats,my_eczema, by.x=c("chromosome", "position"), by.y=c("CHR", "POS"))
#Merge the GWAS and eQTL file by RSID
merged <- merge(summary_stats,my_eczema, by.x="rsid", by.y="RSID")

#Check if the same alleles in both
to_drop <- which((merged$ref != merged$EFFECT_ALLELE & merged$ref != merged$NON_EFFECT_ALLELE) | (merged$alt != merged$EFFECT_ALLELE & merged$alt != merged$NON_EFFECT_ALLELE))
#Drop those rows
if (length(to_drop) > 0){
  merged <- merged[-to_drop,]
}
#Make sure that alleles are harmonised (but not needed for coloc)
#In eQTL Catalog, tte ALT allele is always the effect allele
effect_diff <- which(merged$EFFECT_ALLELE != merged$ALT) # The position of SNPs where effect alleles are different
#Where the alleles are different, flip the beta in eQTL catalog data
merged$beta[effect_diff] <- merged$beta[effect_diff]*(-1)

in_ld <- paste(my_rsid, "_ld.ld", sep="")
in_snps <- paste(my_rsid, "_ld.snplist", sep="")
#Read in LD matrix
my_matrix <- read.delim(in_ld, stringsAsFactors = F, header=F, sep="")
#Read in SNPs present in the matrix
my_snp_list <- read.delim(in_snps, stringsAsFactors = F, header=F, sep="")
colnames(my_matrix) <- my_snp_list$V1
rownames(my_matrix) <- my_snp_list$V1
my_matrix <- data.matrix(my_matrix, rownames.force = NA)
#Subset columns and rows to SNPs present both in eQTL and GWAS.
rsids_to_keep <- intersect(merged$rsid, my_snp_list$V1)
my_matrix <- my_matrix[,rsids_to_keep]
my_matrix <- my_matrix[rsids_to_keep,]
#Sort the merged dataframe in the same order as matrix
merged <- merged[match(colnames(my_matrix), merged$rsid),]
merged <- merged[complete.cases(merged), ]


#Run coloc
tryCatch({
    coloc_res <-  run_coloc(merged, my_matrix)
    },
    error = function(e) {
      #what should be done in case of exeption?
    missing_runs[nrow(missing_runs) + 1,] = c(my_rsid, my_gene, my_hugo_name, row$study, row$tissue_label, row$condition_label, row$quant_method)
    }
  )

#Add gene ENGS, HUGO name, dataset and rsid
coloc_res <- coloc_res %>% tibble::add_column("hugo_name" = my_hugo_name, "ensembl_id" = my_gene, "rsid" = my_rsid,
                                 "study" = row$study, "tissue" = row$qtl_group, 
                                 "condition" = row$condition_label, "feature" = row$quant_method)

#Add to the data frame with all results
coloc_all <- bind_rows(coloc_all, coloc_res)

#Save df to file
output_file <- paste(my_rsid, "susie__coloc.txt", sep="")
output_file2 <- paste(my_rsid, "susie_missing_runs.txt", sep="")
write.table(coloc_all, file=output_file, quote=F, sep="\t", row.names=F)
write.table(missing_runs, file=output_file2, quote=F, sep="\t", row.names=F)
}}