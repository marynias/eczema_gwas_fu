library(gplots)
library(RColorBrewer)
library(tools)

args = commandArgs(trailingOnly=TRUE)

files <- list.files(pattern="*chr*")
txt_files_df <- lapply(files, function(x) {read.table(file = x, header = T, stringsAsFactors=FALSE, row.names=1)})
combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 

df3 <- t(combined_df)

#Create basic palette
mypalette<-brewer.pal(9, "RdPu")
#Ramp up the color palette
cols <- colorRampPalette (mypalette) (20)
par(mar=c(7,2,2,2)+0.1) 
output <- paste(args[1], ".pdf", sep="")
output2 <- paste(args[1], ".txt", sep="")
pdf(output, as.integer(args[2]), as.integer(args[3]))
write.table(df3, output2, quote=F, sep="\t", col.names=NA)
heatmap.2(as.matrix(df3),dendrogram="col", na.color="gray",
          notecol="black",col=cols,scale="none",key=TRUE, keysize=0.75, margins=c(20,15),
          density.info="none", trace="none", cexRow=1, cexCol=1, 
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="single"))
dev.off()

files_paintor <- list.files(pattern="^paintor_chr*")
print(files_paintor)
txt_files_df_paintor <- lapply(files_paintor, function(x) {read.table(file = x, header = T, stringsAsFactors=FALSE, row.names=1)})
df_paintor <- do.call("rbind", lapply(txt_files_df_paintor, as.data.frame)) 
output_paintor <- paste("paintor", "_", args[1], ".txt", sep="")
write.table(df_paintor, output_paintor, quote=F, sep="\t", col.names=NA)


files_finemap <- list.files(pattern="^finemap_chr*")
txt_files_df_finemap <- lapply(files_finemap, function(x) {read.table(file = x, header = T, stringsAsFactors=FALSE, row.names=1)})
df_finemap <- do.call("rbind", lapply(txt_files_df_finemap, as.data.frame)) 
output_finemap <- paste("finemap", "_", args[1], ".txt", sep="")
write.table(df_finemap, output_finemap, quote=F, sep="\t", col.names=NA)

files_jam <- list.files(pattern="^jam_chr*")
txt_files_df_jam <- lapply(files_jam, function(x) {read.table(file = x, header = T, stringsAsFactors=FALSE, row.names=1)})
df_jam <- do.call("rbind", lapply(txt_files_df_jam, as.data.frame)) 
output_jam <- paste("jam", "_", args[1], ".txt", sep="")
write.table(df_jam, output_jam, quote=F, sep="\t", col.names=NA)