library(GenomicRanges)
library(fitCons.UCSC.hg19)
ls("package:fitCons.UCSC.hg19")
fitcons <- fitCons.UCSC.hg19
citation(fitcons)
queries <- read.delim("/panfs/panasas01/sscm/qh18484/analysis/annotation/data_manipulation/interval_r2_0.2_1k_individual.granges", header=FALSE)

res <- lapply(queries$V1, 
  function(x) 
  {
  as.data.frame(gscores(fitcons, GRanges(x)))
  })

all_results <- Reduce(function(x, y) merge(x, y, all=TRUE), res)
annotation <- read.delim("/panfs/panasas01/sscm/qh18484/analysis/annotation/data_manipulation/interval_r2_0.2_1k_individual_padded.bed", header=FALSE)

#Merge rsid annotation with results
colnames(annotation) <- c("seqnames", "start", "end", "rsid")
combined <- merge(all_results, annotation, by=c("seqnames", "start", "end"))
out <- write.table(combined, "interval_r2_0.2_1k_individual.fitCons", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
