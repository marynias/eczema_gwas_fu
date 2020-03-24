library(gplots)
library(RColorBrewer)
library(tools)

args = commandArgs(trailingOnly=TRUE)

files <- list.files()
txt_files_df <- lapply(files, function(x) {read.table(file = x, header = T, stringsAsFactors=FALSE, row.names=1)})
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 

df3 <- t(combined_df)
df3[is.na(df3)] <- 0
#Create basic palette
mypalette<-brewer.pal(9, "RdPu")
#Ramp up the color palette
cols <- colorRampPalette (mypalette) (20)
par(mar=c(7,2,2,2)+0.1) 
output <- args[1]
pdf(output, as.integer(args[2]), as.integer(args[3]))
heatmap.2(as.matrix(df3),dendrogram="col", na.color="gray",
          notecol="black",col=cols,scale="none",key=TRUE, keysize=0.75, margins=c(20,15),
          density.info="none", trace="none", cexRow=1, cexCol=1, 
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="single"))
dev.off()