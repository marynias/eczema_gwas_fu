library(ieugwasr)
library(coloc)
library(dplyr)
library(tools)
library(Hmisc)
library(hablar)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least 3 arguments must be supplied.n", call.=FALSE)}

my_interval <- args[1]
gene_list <- args[2]
gwas_name <- args[3]
#Load required functions
source("/mnt/storage/home/qh18484/bin/eczema_gwas_fu/new_gwas/colocalisation/eqtl_catalogue/ieugwasr_helper_functions.R")

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

missing_runs <- data.frame(rsid = character(),
                           gene_ensembl_id = character(),
                           hugo_name = character(),
                           study = character(),
                           tissue_label = character(),
                           condition_label = character(),
                           quant_method = character(),
                           stringsAsFactors=FALSE
                           
)

coloc_all <-  dplyr::as_tibble(coloc_all)


#Read in the file with target SNPs and interval ranges in GRCh37.
my_ranges <- read.table(my_interval, stringsAsFactors = F, header=F)
colnames(my_ranges) <- c("chrom", "start", "end", "rsid")
my_rsids <- unique(my_ranges$rsid)

#Read in dataframe with hugo symbol and ID mapping to gene name
hugo_gene_ids <- read.delim("hugo_gene_ids.txt", stringsAsFactors = F, header=T, sep="\t")
hugo_gene_ids$Approved.name <- Hmisc::capitalize(hugo_gene_ids$Approved.name)
hugo_gene_ids$HGNC.ID <- sub("HGNC:", "", hugo_gene_ids$HGNC.ID)

#Read in dataframe with UniProt mapping to Hugo symbol and ID.
uniprot_hugo <- read.delim("UniProt_IDs_HGNC.txt", stringsAsFactors = F, header=T, sep="\t")

#Read in dataframe with UniProt mapping to UniProt long names
uniprot_long <- read.delim("UniProt_IDs_long.txt", stringsAsFactors = F, header=T, sep="\t")

#Read in file with the list of transcripts within 1 Mbp of each index SNP.
my_genes <- read.table(gene_list, stringsAsFactors = F, header=F)
colnames(my_genes) <- c("chrom", "start", "end", "rsid", "gene_chrom", "gene_start", "gene_end", "strand", "type", "Ensembl_transcript_id", "Ensembl_gene_id", "hugo_gene_name", "hugo_gene_id", "file")

#IEU GWAS catalog
c <- gwasinfo()
#Filter to only eQTLgen
eqtlgen <- c[c$author=="Vosa U",]
#Filter to only Sun pQTL
sun <- c[c$author=="Sun BB",]
#sun$trait <- gsub("-", "", sun$trait)

merged_genes <- merge(my_genes, hugo_gene_ids, by.x="hugo_gene_id", by.y="HGNC.ID")
#Merge Sun with UniProt IDs.
sun_uniprot <- merge(sun, uniprot_long, by.x="trait", by.y="Target.fullname")
#Merge Sun UniProt IDs with HUGO IDs.
sun_hugo <- merge(sun_uniprot, uniprot_hugo, by.x="UniProt", by.y="UniProt.ID.supplied.by.UniProt.")
#Now merge our genes with Sun traits
merged_sun_0 <- merge(merged_genes, sun_hugo, by.x="Approved.symbol")
#Now merge our genes with eQTLgen traits
merged_eqtlgen_0 <- merge(my_genes, eqtlgen, by.x="Ensembl_gene_id", by.y="trait")

for (my_rsid in my_rsids) {

#Lookup my LD interval range. Need to use a file with GRCh38 coordinates.
my_chrom = my_ranges[my_ranges$rsid == my_rsid,]$chrom
my_chrom = sub("chr", "", my_chrom)
my_start = my_ranges[my_ranges$rsid == my_rsid,]$start
my_end= my_ranges[my_ranges$rsid == my_rsid,]$end
coordinates <- paste(my_chrom, ":", my_start, "-", my_end, sep="")

#Load file with my eczema summary stats for a given variant.
eczema_file = paste("GWAS_intervals/", my_rsid, "_", gwas_name, ".gwas",sep="")
my_eczema <- read.table(eczema_file, stringsAsFactors = F, header=T)
my_eczema <- my_eczema %>% dplyr::as_tibble() %>% dplyr::rename(MAF = maf) %>%
  hablar::convert(num(CHR, POS, N, MAF, BETA, SE, Z_SCORE, PVAL))

#Find proteins which overlap between Sun and our Dataset.
#Find genes of interest for current locus for Sun
merged_sun = merged_sun_0[merged_sun_0$rsid == my_rsid,]
#Find genes of interest for current locus for eQTLgen
merged_eqtlgen = merged_eqtlgen_0[merged_eqtlgen_0$rsid == my_rsid,]

if (dim(merged_sun)[1] > 0) {
process_coloc(merged_sun, my_eczema, "Sun2018", "protein", coordinates, gwas_name)}
if (dim(merged_eqtlgen)[1] > 0) {
process_coloc(merged_eqtlgen, my_eczema, "eQTLgen", "ge", coordinates, gwas_name)}
}
