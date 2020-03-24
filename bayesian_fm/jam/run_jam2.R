# By default 1 million iterations are run

library(R2BGLiMS) # Load package
library(tools)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("At least three arguments must be supplied (2 input files and Priority Pruner results).n", call.=FALSE)}

my_betas <- args[1]
my_dosage <- args[2]
pp_results <- args[3]


my_betas_df <- read.table(my_betas, header=T, stringsAsFactors=FALSE)
my_dosage_df <- read.table(my_dosage, header=T, row.names=1, stringsAsFactors=FALSE)
pp_results_df <- read.table(pp_results, header=T, stringsAsFactors=FALSE)
#Filter only SNPs with selected == 1 (ie. not pruned)
pp_results_sel <- subset(pp_results_df, selected == "1")
my_betas_df <- my_betas_df[my_betas_df$rsid %in% pp_results_sel$name, ]
my_dosage_df <- subset(my_dosage_df, rownames(my_dosage_df) %in% pp_results_sel$name)

my_dosage_df <- t(my_dosage_df)

#Estimate missing values as the mean value of dosage across all the indvs for the snp, as suggested by Paul.
for(i in 1:ncol(my_dosage_df)){
  my_dosage_df[is.na(my_dosage_df[,i]), i] <- mean(my_dosage_df[,i], na.rm = TRUE)
}

my_dosage_matrix <- as.matrix(my_dosage_df)

#Number of samples in eczema GWAS
n=103066
snp.sds = apply(my_dosage_matrix,MAR=2,sd)
#Proportion of cases in eczema GWAS
p.cases=0.2245

# Standardised effects (z-scores)
my_betas_df$zscores <- my_betas_df$beta / (my_betas_df$se*sqrt(n))
# Divide by SNP SDs to get allelic effects
my_betas_df$zscores <- my_betas_df$zscores / snp.sds # Divide standardised linear effects by SNP standard deviations
# Adjust for fraction of cases
my_betas_df$zscores <- my_betas_df$zscores*sqrt(p.cases*(1-p.cases)) # Multiply by trait SD for effect on trait scale


my_betas_vector <- structure(as.numeric(my_betas_df$zscore), names = as.character(my_betas_df$rsid))

extra_args <- list(
  "GaussianResidualVarianceInvGammaPrior_a" = 2,
  "GaussianResidualVarianceInvGammaPrior_b" = p.cases*(1-p.cases)
)


jam.results <- JAM(
  marginal.betas=my_betas_vector,
  X.ref=my_dosage_matrix,
  model.space.prior = list("a"=1, "b"=length(names(my_betas_vector)), "Variables"=names(my_betas_vector)),
  n.mil.iter = 10,
  trait.variance = p.cases*(1-p.cases), # Optional: Invokes a new refinement to JAM
  extra.arguments = extra_args,
  n=n)

output <- paste(file_path_sans_ext(my_betas), "_manhattan.pdf", sep="")
pdf(output)
ManhattanPlot(jam.results, plot.title=file_path_sans_ext(my_betas))
dev.off()

output2 <- paste(file_path_sans_ext(my_betas), "_full_table.txt", sep="")
my_results <- jam.results@posterior.summary.table
write.table(my_results, file=output2, sep="\t",quote=FALSE, col.names=NA, row.names=TRUE)

output3 <- paste(file_path_sans_ext(my_betas), "_pretty_table.txt", sep="")
my_pretty <- PrettyResultsTable(jam.results, round.digits.betas = 3,
                   round.digits.postprob = 3, round.digits.bf = 3, normalised.sds = NULL)
write.table(my_pretty, file=output3, sep="\t",quote=FALSE, col.names=NA, row.names=TRUE)

TopModels(jam.results) 

output4 <- paste(file_path_sans_ext(my_betas), "_credible_set.txt", sep="")
cs <- CredibleSet(jam.results, credible.percentile.threshold = 0.95)
write.table(cs, file=output4, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)

output5 <- paste(file_path_sans_ext(my_betas), "_region_bf.txt", sep="")
mbf <- ModelSizeBayesFactors(jam.results)
write.table(mbf, file=output5, sep="\t",quote=FALSE, col.names=NA, row.names=TRUE)