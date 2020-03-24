library("dmGWAS")
geneweight <- read.delim("results_euro_pval_1k_vegas.dmgwas", header=F)
colnames(geneweight) <- c("gene", "weight")
network <- read.delim("9606.protein.links.full.v10.5.conf.dmgwas", header=F)
colnames(network) <- c("interactorA", "interactorB")
#First, cases=lesional skin
expr1 <- read.delim("Suarez-Farinas2011/all_lesion.expr", header=T)
colnames(expr1)[1] <-"Gene_symbol"
#Control=normal
expr2 <- read.delim("Suarez-Farinas2011/all_normal.expr", header=T)
colnames(expr2)[1] <-"Gene_symbol"

#Subset geneweights to genes which are in the network.
geneweight <- geneweight[geneweight$gene %in% network$interactorA | geneweight$gene %in% network$interactorB, ]
#Subset geneweights to genes which have expression data.
#geneweight <- geneweight[geneweight$gene %in% expr1$Gene_symbol, ]

#Subset network to genes which have expression data.
#network <- network[network$interactorB %in% expr1$Gene_symbol, ]
#network <- network[network$interactorA %in% expr1$Gene_symbol, ]

#Substitute all the 1 values with 0.9999999999999999.
geneweight[geneweight$weight==1,c("weight")] <- 0.9999999999999999

#res.list <- dms(network, geneweight, expr1, expr2, r=0.1, lambda="default")
#selected <- chooseModule(res.list, top=0.01, plot=T)

res.list <- dms(network, geneweight, expr1=NULL, expr2=NULL, d=1, r=0.1)
selected <- chooseModule(res.list, top=0.01, plot=T)


#Now same analysis but contrasting non-lesional AD patients' skin with normal skin
#expr1 <- read.delim("Suarez-Farinas2011/all_nonlesion.expr", header=T)
#colnames(expr1)[1] <-"Gene_symbol"

#res.list2 <- dms(network, geneweight, expr1, expr2, r=0.1, lambda="default")
#selected2 <- chooseModule(res.list2, top=0.01, plot=T)

#save(res.list, res.list2, selected, selected2, file = "my_dmgwas.RData")
save(res.list, selected, file = "my_dmgwas.RData")