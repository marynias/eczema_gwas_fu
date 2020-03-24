library(gpart)
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (2 input files and Priority Pruner results).n", call.=FALSE)}

my_snp <- args[1]
my_dosage <- args[2]
heuristic <- args[3]

geno <- read.table(my_dosage, header=T, stringsAsFactors=FALSE, row.names=1) 
geno <- t(geno)
SNPinfo <- read.table(my_snp, header=T, stringsAsFactors=FALSE)
SNPinfo <- cbind(chrN = "Z" , SNPinfo)

#CLQD partitioning the SNPs into subgroups such that each subgroup contains highly correlated SNPs.
#CLQres = CLQD(geno, SNPinfo, CLQmode = 'Density')

#Big_LD returns the estimation of LD block regions of given data.
BigLDres = BigLD(geno, SNPinfo, LD = "Dprime", hrstType = heuristic)
output_file <-  paste(file_path_sans_ext(my_dosage), ".gpart", sep="")
write.table(BigLDres, file=output_file, sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)

