library(gpart)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (2 input files and Priority Pruner results).n", call.=FALSE)}

my_snp <- args[1]
my_dosage <- args[2]
results <- args[3]
my_chrom <- args[4]
snp_id <- args[5]
my_position <- args[6]
my_interval <- args[7]
my_dataset <- args[8]


geno <- read.table(my_dosage, header=T, stringsAsFactors=FALSE, row.names=1) 
geno <- t(geno)
SNPinfo <- read.table(my_snp, header=T, stringsAsFactors=FALSE)
SNPinfo <- cbind(chrN = "Z" , SNPinfo)
my_results <- read.table(results, header=T, stringsAsFactors=FALSE) 

my_start <- as.integer(my_position) - as.integer(my_interval)
my_end <- as.integer(my_position) + as.integer(my_interval)

my_filename <- paste("chr", my_chrom, snp_id, my_position, my_interval, my_dataset, "ldplot", sep=".")

LDblockHeatmap(geno = geno, SNPinfo = SNPinfo, blockresult=my_results, startbp = my_start, endbp = my_end, 
               filename = my_filename, res=300)
