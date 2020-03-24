library("rtracklayer")
my_human <- import.gff("Homo_sapiens.GRCh37.87.gff3.gz")
my_human_genes <- my_human$type %in% c("gene", "lincRNA_gene", "rRNA_gene", "snRNA_gene", "snoRNA_gene", "miRNA_gene", "mt_gene", "processed_transcript", "processed_pseudogene", "RNA", "V_gene_segment", "VD_gene_segment", "rRNA", "J_gene_segment", "pseudogene", "miRNA", "snoRNA", "C_gene_segment")
gene_data <- my_human[my_human_genes,]
#Do we have a gene ID?
my_human_genes2 <- !is.na(gene_data$gene_id)
gene_data <- gene_data[my_human_genes2,]
my_range <- import.bed("paternoster_2015_index_snps_sorted_3Mbp_nodups.bed")
hts_within <- findOverlaps(gene_data, my_range, type='within', select="all", ignore.strand=TRUE)
my_range$gene_count <-countOverlaps(my_range,gene_data, ignore.strand=TRUE)

#Counts of number of genes per index locus 3Mbp interval.
df2 <- data.frame(chr=seqnames(my_range),
                  starts=start(my_range),
                  ends=end(my_range),
                  snp_name=my_range$name,
                  gene_count=my_range$gene_count)
write.table(df2, file="paternoster_2015_index_snps_sorted_3Mbp_nodup.gene_counts", quote=F, sep="\t", col.names=T, row.names=F)

range_results <- my_range[subjectHits(hts_within)]
gene_results <- gene_data[queryHits(hts_within)]
range_results$ensembl <- sub("^gene:", "", gene_results$ID)
range_results$gene_name <- gene_results$Name

#Each index SNP matched to all the genes within the 3Mbp interval around it.
df <- data.frame(chr=seqnames(range_results),
                 starts=start(range_results),
                 ends=end(range_results),
                 snp_name=range_results$name,
                 ensembl=range_results$ensembl,
                 gene_name=range_results$gene_name
                 )
write.table(df, file="paternoster_2015_index_snps_sorted_3Mbp_nodup.ensembl.gene_matches", quote=F, sep="\t", col.names=T, row.names=F)

