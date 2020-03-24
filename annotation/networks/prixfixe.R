library(PrixFixe, quietly = TRUE)
gwas <- read.table("non_redundant_loci.txt", header=F, stringsAsFactors=FALSE)
ld_regions_query_result <- getLDRegionsFromSNPs.WS(gwas$V1)
candidate_gene_table <- getGeneRegionsTable(ld_regions_query_result$ld_regions)
candidate_edges <- getCandidateEdges.WS(candidate_gene_table$id)
candidate_gene_sets <- tapply(candidate_gene_table$id, candidate_gene_table$region_name, I)
prix_fixe <- GPF$PF$new(candidate_gene_sets, candidate_edges, VERBOSE = TRUE)
GA_run_result <- GPF$GA$run(prix_fixe, VERBOSE = TRUE)
score_table <- GPF$GA$scoreVertices(GA_run_result, VERBOSE = TRUE)
candidate_gene_table$score <- score_table$scaled_score[match(candidate_gene_table$id, score_table$vertex_name)]
candidate_gene_table$score[is.na(candidate_gene_table$score)] <- 0
library(plyr)
candidate_gene_table <- plyr::arrange(candidate_gene_table, region_name, -score)
candidate_gene_table[c("id", "symbol", "region_name", "score")]
loci_mapping <- read.table("loci_mapping.txt", header=T, stringsAsFactors=TRUE)
id <- loci_mapping$rsid[match(candidate_gene_table$region_name, loci_mapping$locus)]
candidate_gene_table$rsid <- id
write.table(candidate_gene_table, file="non_redundant_loci.prixfixe.out", quote=F, sep="\t", col.names=T)

