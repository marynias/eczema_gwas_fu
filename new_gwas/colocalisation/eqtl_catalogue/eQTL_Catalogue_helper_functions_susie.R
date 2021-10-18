library("dplyr")
library("coloc")

import_eQTLCatalogue_no_rsid <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics
  tbx <- TabixFile(ftp_path)
  res <- scanTabix(tbx, param=region)
  fetch_table  <- res %>% dplyr::as_tibble() %>% dplyr::rename(value = 1) %>% tidyr::separate("value", column_names, "\t")
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats <- dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::mutate_at(vars(chromosome, position, maf, pvalue, beta, se, ac, an), list(as.numeric)) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
  return(summary_stats)
}

import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics
  tbx <- TabixFile(ftp_path)
  res <- scanTabix(tbx, param=region)
  res <- scanTabix(ftp_path, param=region)
  fetch_table  <- res %>% dplyr::as_tibble() %>% dplyr::rename(value = 1) %>% tidyr::separate("value", column_names, "\t")
  summary_stats <- dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::mutate_at(vars(chromosome, position, maf, pvalue, beta, se, ac, an), list(as.numeric)) %>%
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
  return(summary_stats)
}

import_eQTLCatalogue_ver_tabix <- function(ftp_path, coordinates, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics
  my_command <- paste("tabix ", ftp_path, " ", coordinates, sep="")
  #Try to get summary stats 10 times before giving up

  fetch_table <-  system(my_command, intern=TRUE, timeout=600) %>% dplyr::as_tibble() %>% dplyr::rename(value = 1) %>% tidyr::separate("value", column_names, "\t")
  summary_stats <- dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::mutate_at(vars(chromosome, position, maf, pvalue, beta, se, ac, an), list(as.numeric)) %>%
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1)  #Multiallelics
#NA
  return(summary_stats)
}


import_eQTLCatalogue_ver_tabix_tx <- function(ftp_path, coordinates, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics
  my_command <- paste("tabix ", ftp_path, " ", coordinates, sep="")
  #Try to get summary stats 10 times before giving up

  fetch_table <-  system(my_command, intern=TRUE, timeout=600) %>% dplyr::as_tibble() %>% dplyr::rename(value = 1) %>% tidyr::separate("value", column_names, "\t")
  summary_stats <- dplyr::filter(fetch_table, molecular_trait_id == selected_gene_id) %>%
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::mutate_at(vars(chromosome, position, maf, pvalue, beta, se, ac, an), list(as.numeric)) %>%
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1)  #Multialllics
 #NA
  return(summary_stats)

}

#SNPs identified by their RSIDs here.
run_coloc <- function(merged, my_matrix){
  eQTL_dataset = list(beta = merged$beta,
                      varbeta =  merged$SE^2,
                      N = (merged$an)[1]/2, # Samples size is allele number (AN) dvided by 2
                      MAF = merged$maf, 
                      type = "quant", 
                      snp = merged$rsid,
                      LD = my_matrix)
  gwas_dataset = list(beta = merged$BETA,
                      varbeta = merged$SE^2, 
                      type = "cc", 
                      snp = merged$rsid,
                      MAF = merged$maf, 
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


