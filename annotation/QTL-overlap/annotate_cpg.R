library(tools)
library(FDb.InfiniumMethylation.hg19)

args = commandArgs(trailingOnly=TRUE)

my_file <- args[1]

cpgs <- read.table(my_file, header=F, stringsAsFactors=FALSE)
header <- c('snp', 'rsid', 'cpg', 'Allele1', 'Allele2', 'Freq1', 'FreqSE', 'Effect', 'StdErr', 'Pvalue', 'Direction', 'HetISq',	'HetChiSq',
'HetDf', 'HetPVal',	'EffectARE', 'StdErrARE', 'PvalueARE', 'tausq',	'StdErrMRE', 'PvalueMRE', 'TotalSampleSize', 'snpchr', 'snppos',
'snptype',	'cpgchr', 'cpgpos',	'cis')
colnames(cpgs) <- header
hm450 <- get450k()
my_query <- unique(cpgs$cpg)

probes <- hm450[my_query]
tss <- getNearestTSS(probes)
transcript <- getNearestTranscript(probes)

#Merge with CpG Ids.
tss$cpg <- my_query
transcript$cpg <- my_query
names(tss)[names(tss) == 'nearestGeneSymbol'] <- 'nearest_tss'
names(transcript)[names(transcript) == 'nearestGeneSymbol'] <- 'nearest_transcript'
tss_slim <- tss[,c("cpg", "nearest_tss")]
transcript_slim <- transcript[,c("cpg", "nearest_transcript")]
annotation <- merge(tss_slim, transcript_slim, by="cpg")
all <- merge(cpgs, annotation, by.x="cpg", by.y="cpg", all.x=TRUE)

#Merge with input.
write.table(all, paste(my_file, "_genes", sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)