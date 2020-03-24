library(tools)
library(GenABEL)

#Read in the covariate file, with disease status and age (gender is the same for all samples)
my_covar <- read.delim("twinsuk_covariates2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
my_covar <- as.data.frame(t(my_covar))

#Read in file mapping ensembl_ids to gene names.
ensembl <- read.delim("ensemble_ids_unique_ensembl_names", sep="\t", stringsAsFactors=FALSE, header=TRUE)# load file

files <- list.files(pattern="*_headers_eczema.rpkm$", recursive=F)
lapply(files, function(x) {
  my_rpkm_df <- read.delim(x, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
  #Subset to individuals with genotypes and reorder to genotype order.
  my_order <-  read.delim("genosample_eczema", sep="\t", stringsAsFactors=FALSE, header=FALSE)
  my_ids <- unlist(my_order, use.names=FALSE)
  #First drop not needed individuals
  my_rpkm_df <- my_rpkm_df[my_ids]
  #Same but for covariates
  my_covar <- my_covar[my_ids]
  my_covar <- t(my_covar)
  #Drop the disease status column
  my_covar <- my_covar[, c(-2)]
  #Write covariate to file.
  write.table(my_covar, file="gemma_eczema_covariate.txt", quote=F, sep="\t", col.names = F, row.names = F)
  for(i in 1:nrow(my_rpkm_df)) 
  {
    my_row <- as.numeric(my_rpkm_df[i,])
    if (length(unique(my_row)) < 3) 
    {
      my_rpkm_df[i,] <- my_row
    }
    else
    {
      my_rpkm_df[i,] <- rntransform(my_row)
      my_gene_id <- rownames(my_rpkm_df[i,])
      my_genes <- ensembl[ensembl$gene==my_gene_id,]
      my_gene = my_genes$gene_names
      transposed_df <- t(my_rpkm_df[i,])
      out <- paste(my_gene, "_", x, ".rnt", sep="")
      write.table(transposed_df, file=out, quote=F, sep="\t", col.names = F, row.names = F)
    }
  }
  
})


my_covar <- read.delim("twinsuk_covariates2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
my_covar <- as.data.frame(t(my_covar))

files <- list.files(pattern="*_headers_noneczema.rpkm$", recursive=F)
lapply(files, function(x) {
  my_rpkm_df <- read.delim(x, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
  #Subset to individuals with genotypes and reorder to genotype order.
  my_order <-  read.delim("genosample_noneczema", sep="\t", stringsAsFactors=FALSE, header=FALSE)
  my_ids <- unlist(my_order, use.names=FALSE)
  #First drop not needed individuals
  my_rpkm_df <- my_rpkm_df[my_ids]
  #Same but for covariates
  my_covar <- my_covar[my_ids]
  my_covar <- t(my_covar)
  #Drop the disease status column
  my_covar <- my_covar[, c(-2)]
  #Write covariate to file.
  write.table(my_covar, file="gemma_noneczema_covariate.txt", quote=F, sep="\t", col.names = F, row.names = F)
  for(i in 1:nrow(my_rpkm_df)) 
  {
    my_row <- as.numeric(my_rpkm_df[i,])
    if (length(unique(my_row)) < 3) 
    {
      my_rpkm_df[i,] <- my_row
    }
    else
    {
      my_rpkm_df[i,] <- rntransform(my_row)
      my_gene_id <- rownames(my_rpkm_df[i,])
      my_genes <- ensembl[ensembl$gene==my_gene_id,]
      my_gene = my_genes$gene_names
      transposed_df <- t(my_rpkm_df[i,])
      out <- paste(my_gene, "_", x, ".rnt", sep="")
      write.table(transposed_df, file=out, quote=F, sep="\t", col.names = F, row.names = F)
    }
  }
  
})

my_covar <- read.delim("twinsuk_covariates2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
my_covar <- as.data.frame(t(my_covar))

#When it comes to all samples, add disease status as phenotype.
files <- list.files(pattern="*_headers_all.rpkm$", recursive=F)
lapply(files, function(x) {
  my_rpkm_df <- read.delim(x, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)# load file
  #Subset to individuals with genotypes and reorder to genotype order.
  my_order <-  read.delim("genosample_all", sep="\t", stringsAsFactors=FALSE, header=FALSE)
  my_ids <- unlist(my_order, use.names=FALSE)
  #First drop not needed individuals
  my_rpkm_df <- my_rpkm_df[my_ids]
  #Same but for covariates
  my_covar <- my_covar[my_ids]
  my_covar <- t(my_covar)
  #Save the disease status column
  my_disease <-  my_covar[, c(-1)]
  #Drop the disease status column
  my_covar <- my_covar[, c(-2)]
  #Write covariate to file.
  write.table(my_covar, file="gemma_all_covariate.txt", quote=F, sep="\t", col.names = F, row.names = F)
  for(i in 1:nrow(my_rpkm_df)) 
  {
    my_row <- as.numeric(my_rpkm_df[i,])
    if (length(unique(my_row)) < 3) 
    {
      my_rpkm_df[i,] <- my_row
    }
    else
    {
      my_rpkm_df[i,] <- rntransform(my_row)
      my_gene_id <- rownames(my_rpkm_df[i,])
      my_genes <- ensembl[ensembl$gene==my_gene_id,]
      my_gene = my_genes$gene_names
      transposed_df <- t(my_rpkm_df[i,])
      #Bind the column with disease status
      transposed_df2 <- cbind(my_disease, transposed_df)
      out <- paste(my_gene, "_", x, ".rnt", sep="")
      write.table(transposed_df2, file=out, quote=F, sep="\t", col.names = F, row.names = F)
    }
  }
})



