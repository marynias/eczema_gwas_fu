#SNPs identified by their RSIDs here.
run_coloc <- function(merged){
  eQTL_dataset = list(beta = merged$beta,
                      varbeta = merged$se^2,
                      N = merged$n, # Samples size is allele number (AN) dvided by 2
                      MAF = merged$MAF, 
                      type = "quant", 
                      snp = merged$rsid)
  gwas_dataset = list(beta = merged$BETA,
                      varbeta = merged$SE^2, 
                      type = "cc", 
                      snp = merged$rsid,
                      MAF = merged$MAF, 
                      N = merged$N,
                      s=0.2245)
  #Recommended prior values. P12 is a good choice as only testing relevant tissues shown to be enriched.
  coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  res_formatted = dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
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
    coloc_res <- run_coloc(merged)
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
  output_file <- paste(my_rsid, "_", my_study, "_ieu_coloc.txt", sep="")
  output_file2 <- paste(my_rsid, "_", my_study, "_ieu_missing_runs.txt", sep="")
  write.table(coloc_all, file=output_file, quote=F, sep="\t", row.names=F)
  write.table(missing_runs, file=output_file2, quote=F, sep="\t", row.names=F)
  
}}



