########## moloc on single locus

#The summary statistics for each dataset must be in a list (e.g. list(gwas, eqtl, mqtl)).
#Must have columns `SNP`, `BETA`, `SE`.

library(knitr)
options(scipen = 1, digits = 2)
## load single locus data (in a list) and bed file
data_single=get(load(file="data_single.rda"))
t1 = knitr::kable(head(data_single[[1]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t2 = knitr::kable(head(data_single[[2]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t3 = knitr::kable(head(data_single[[3]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
cat(c('<table><tr valign="top"><td>', t1, '</td><td>', t2, '</td><td>', t3, '</td><tr></table>'),
    sep = '')

options(scipen = 1, digits = 2)
library(moloc)
moloc <- moloc_test(data_single, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07))
# Posteriors
print(moloc[[1]])
# Number of SNPs analyzed
print(moloc[[2]])
# Posterior of the most likely SNP co-localizing with another trait
print(moloc[[3]])


########## moloc Genome-wide/multiple loci analysis

#The summary statistics for each dataset in a list (e.g. list(gwas, eqtl, mqtl)), same as for the single locus, but now have multiple ProbeIDs;
#Must have columns `SNP`, `BETA`, `SE`; 
#If the regions are defined based on ProbeID, these must match the ProbeID in the file.
#The bed file: specifies the region to use.
#Must contain the columns `ProbeID`, `CHR`, `START`, `STOP` (other columns are ignored)

library(knitr)
options(scipen = 1, digits = 2)
## load single locus data (in a list) and bed file
library(moloc)
data_genome=get(load(file="data_genome.rda"))
t1 = knitr::kable(head(data_genome[[1]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t2 = knitr::kable(head(data_genome[[2]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
t3 = knitr::kable(head(data_genome[[3]], 5),  format='html', table.attr='cellpadding="1"', output = FALSE)
cat(c('<table><tr valign="top"><td>', t1, '</td><td>', t2, '</td><td>', t3, '</td><tr></table>'),
    sep = '')
knitr::kable(bed)

options(scipen = 1, digits = 2)
library(moloc)
library(foreach)
library(doParallel)
moloc_genome <- coloc.genome(data_genome, bed, prefix = "pref", save.SNP.info=FALSE, cores=20, have_alleles=TRUE, bychrpos=TRUE, prior_var="default", priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 50, takelog = FALSE, write=TRUE, outfolder = "test", forcePVAL = FALSE)
# Posteriors
print(moloc_genome[[1]])
# Bed file of loci analyzed
print(moloc_genome[[2]])