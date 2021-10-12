library(tools)
library(dplyr)
library(ghql)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("At least 1 argument must be supplied.n", call.=FALSE)}

input_stats <- args[1]
gwas_name <- args[2]


query_v2g <- function(my_variant){
  con <- GraphqlClient$new("http://genetics-api.opentargets.io/graphql")
  qry <- Query$new()
  query_string <- '{
  genesForVariant(variantId: "placeholder") {
    overallScore
    gene {
      id
      symbol
    }
  }
}'
  #Substitute placeholder with our variant of interest.
  query_string <- gsub("placeholder", my_variant, query_string)

  qry$query('variant', query_string)
  res <- con$exec(qry$queries$variant)
  final_table <- jsonlite::fromJSON(res)$data$genesForVariant
  final_table <- as.data.frame(final_table)
  final_table <- cbind(final_table$overallScore, final_table$gene)
  final_table$variant <- my_variant
  #If no result for variant returned
  if (length(final_table) == 1) {return (NULL)}
  return(final_table)
}


query_variant <- function(my_variant) {
con <- GraphqlClient$new("http://genetics-api.opentargets.io/graphql")
qry <- Query$new()
query_string <-'{
  search(queryString:"placeholder"){
    totalVariants
    variants {
      id
    }
  }
}'
query_string <- gsub("placeholder", my_variant, query_string)
qry$query('variant', query_string)
res <- con$exec(qry$queries$variant)
final_table <- jsonlite::fromJSON(res)$data$search
if (final_table$totalVariants == 0) {return (NULL)}
final_table <- as.data.frame(final_table)
final_table$variant <- my_variant
return(final_table)
}


df <- read.delim(input_stats, header=F, stringsAsFactors=F,row.names=NULL, sep="\t")
colnames(df) <- c("CHR", "POS", "RSID")

gwas_variants <- lapply(unique(df$RSID), query_variant)
gwas_input <- do.call(rbind,gwas_variants)
gwas_v2g_results <- lapply(gwas_input$id, query_v2g)
gwas_v2g <- do.call(rbind,gwas_v2g_results)

#Bind to get rsids
gwas_v2g <- merge(gwas_v2g, gwas_input, by.x="variant", by.y="id", all.x=T)
gwas_v2g <- gwas_v2g %>% rename(rsid = variant.y)

output <- paste0(gwas_name, "_open_targets.txt")
write.table(gwas_v2g, output, sep="\t", quote=F, row.names=F, col.names=T)