#SNPs identified by their RSIDs here.
run_coloc <- function(merged, my_matrix){
  eQTL_dataset = list(beta = merged$beta,
                      varbeta = merged$se^2,
                      N = merged$n, # Samples size is allele number (AN) divided by 2
                      MAF = merged$MAF, 
                      type = "quant", 
                      snp = merged$rsid,
                      LD = my_matrix)
  gwas_dataset = list(beta = merged$BETA,
                      varbeta = merged$SE^2, 
                      type = "cc", 
                      snp = merged$rsid,
                      MAF = merged$MAF, 
                      N = merged$N,
                      s=0.2245,
                      LD = my_matrix)
  #Recommended prior values. P12 is a good choice as only testing relevant tissues shown to be enriched.
  tryCatch({
  S3=runsusie(eQTL_dataset,nref=503,check_R = FALSE)
  S4=runsusie(gwas_dataset,nref=503,check_R = FALSE)})
  susie.res=coloc.susie(S3,S4)
  res <- susie.res$summary
  #Keep only row with highest PP of H4
  row_to_keep = which.max(res$PP.H4.abf)
  res <- as_tibble(res[row_to_keep,])
  #Drop unncessary columns
  res_formatted <- res[-c(2,3,9,10),]
  res_formatted <- res_formatted %>% select (-c(hit1, hit2, idx1, idx2))
  return(res_formatted)
}

process_coloc <- function(merged_sun, my_eczema, my_study, my_feature, coordinates) {
for (a in 1:nrow(merged_sun)) {
  row <- merged_sun[a,] 
  my_gwas_id <- row$id
  my_gene <- row$Ensembl_gene_id
  my_hugo_name <- row$hugo_gene_name
  my_tissue <- "whole blood"
  my_condition <- "NA"
  
  #Get summary stats
  tryCatch({
    summary_stats <- ieugwasr::associations(coordinates, my_gwas_id)
    merged <- merge(summary_stats,my_eczema, by.x="rsid", by.y="RSID")
  },
  error = function(e) {
    #what should be done in case of exeption?
    missing_runs[nrow(missing_runs) + 1,] = c(my_rsid, my_gene, my_hugo_name, my_study, my_tissue, my_condition, my_feature)
  }
  )
  
  #Run coloc
  tryCatch({
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
    coloc_res <- run_coloc(merged, my_matrix)
    coloc_res <- coloc_res %>% tibble::add_column("hugo_name" = my_hugo_name, "ensembl_id" = my_gene, "rsid" = my_rsid,
                                                  "study" = my_study, "tissue" = my_tissue, 
                                                  "condition" = my_condition, "feature" = my_feature)
    
  },
  error = function(e) {
    #what should be done in case of exeption?
    missing_runs[nrow(missing_runs) + 1,] = c(my_rsid, my_gene, my_hugo_name, my_study, my_tissue, my_condition, my_feature)
  }
 )
  

  #Add to the data frame with all results
  coloc_all <- bind_rows(coloc_all, coloc_res)
  
  #Save df to file
  output_file <- paste(my_rsid, "_", my_study, "_susie_ieu_coloc.txt", sep="")
  output_file2 <- paste(my_rsid, "_", my_study, "_susie_ieu_missing_runs.txt", sep="")
  write.table(coloc_all, file=output_file, quote=F, sep="\t", row.names=F)
  write.table(missing_runs, file=output_file2, quote=F, sep="\t", row.names=F)
  
}}



