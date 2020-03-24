##########################################################
##### Assoc plot #####
##########################################################

#' assoc_plot3
#'
#' assoc_plot3 plots a scatter graph of associations (e.g. log10 p-values)
#' @param data data.frame with markername (marker), chromosome (chr), position (pos) and either z-statistics (z) or probabilities (prob) columns
#' @param corr correlation matrix between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param ylab the y-axis label
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param type the type of the plot either log10p or probabilities
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
assoc_plot3 <- function(data, corr=NULL, corr.top=NULL, ylab=NULL, title=NULL, subtitle=NULL, type="log10p", x.min=NULL, x.max=NULL, top.marker=NULL, legend=TRUE){
  
  # Error messages
  if(!(type=="log10p" | type=="prob")) stop("the type of plot has to be either log10p or prob")
  if(type=="log10p"){
    if(any(names(data)!=c("marker", "chr", "pos", "z"))) stop("dataset needs to include marker, chr, pos and z columns in that order")
  }else{
    if(any(names(data)!=c("marker", "chr", "pos", "prob"))) stop("dataset needs to include marker, chr, pos and prob columns in that order")
  }
  if(!is.null(corr)){if(ncol(corr)!=nrow(data) | nrow(corr)!=nrow(data)) stop("corr has to have the same dimensions as the number of rows in the markers dataset")}
  # if(any(rownames(corr)!=data$marker)) stop("corr has to have the same markers in the same order as the dataset")
  if(length(unique(data$chr))>1) stop("there should only be markers from one chromosome in the markers dataset") 
  if(!(data$chr[1] %in% 1:22)) stop("the plotting tool is only for autosomal chromosomes") 
  if(any(is.na(data))) stop("there are missing values in the dataset") 
  if(class(data$pos)!="integer") stop("the pos variable has to be an integer")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  if(is.null(corr) & !is.null(corr.top)){if(length(corr.top)!=nrow(data)) stop("corr.top has to have the same length as the number of rows in the markers dataset")}
  if(!is.null(top.marker) & length(which(top.marker==data$marker))==0) stop("top.marker is not contained in the markers dataset")
  if(!is.null(top.marker) & length(which(top.marker==data$marker))>1) stop("top.marker maps to multiple markers in the markers dataset")
  
  # Dataset
  if(type=="log10p"){
    mlog10p <- -(log(2) + pnorm(-abs(data$z), log.p=T))/log(10)
    mlog10p[mlog10p>1000] <- 1000
    data$stats <- mlog10p
  }else{data$stats <- data$prob}
  data <- data[,c("marker", "chr", "pos", "stats")]
  data$marker <- as.character(data$marker)
  chr <- as.integer(data$chr[1])
  if(is.null(x.min)){x.min <- min(as.integer(data$pos))}
  if(is.null(x.max)){x.max <- max(as.integer(data$pos))}
  if((x.max - x.min)>10000000) stop("the plotting tool can plot a maximum of 10MB")
  
  # Genes
  gene.region <- genes[genes$chr==chr & !(genes$end<x.min) & !(genes$start>x.max),]
  gene.region$start[gene.region$start<x.min] <- x.min
  gene.region$end[gene.region$end>x.max] <- x.max
  gene.region <- gene.region[with(gene.region, order(start)), ]
  ngenes <- nrow(gene.region)
  
  # Max and min
  x.min <- x.min - 0.02*(x.max - x.min)
  x.max <- x.max + 0.02*(x.max - x.min)
  
  # Correlation matrix
  if(is.null(corr) & is.null(corr.top)){r2_legend <- FALSE; corr <- matrix(NA, nrow=nrow(data), ncol=nrow(data))}
  
  # Recombination plot
  recombination.plot <- plot_recombination_rate(chr, x.min, x.max)
  
  # Gene plot
  if(ngenes==0){gene.plot <- plot_gene_zero(chr, x.min, x.max)}
  if(ngenes>0 & ngenes<=5){gene.plot <- plot_gene_two(gene.region, chr, x.min, x.max)}
  if(ngenes>5 & ngenes<=10){gene.plot <- plot_gene_five(gene.region, chr, x.min, x.max)}
  if(ngenes>10 & ngenes<=25){gene.plot <- plot_gene_ten(gene.region, chr, x.min, x.max)}
  if(ngenes>25){gene.plot <- plot_gene_fifteen(gene.region, chr, x.min, x.max)}
  
  # Marker plot
  data$chr <- as.integer(data$chr)
  data$pos <- as.integer(data$pos)
  if(type=="log10p"){ylab <- expression("-log"["10"]*paste("(",italic("p"),")"))}else{if(is.null(ylab)){ylab <- "Probability"}}  
  marker.plot <- plot_assoc3(data, corr, corr.top, x.min, x.max, top.marker, ylab, type)
  
  # Combined plot
  legend <- FALSE
  combined.plot <- plot_assoc_combined3(recombination.plot, gene.plot, marker.plot, title, subtitle, ngenes, legend)
  
  return(combined.plot)
}
##########################################################
##### Combined plot #####
##########################################################

#' plot_assoc_combined3
#'
#' plot_assoc_combined combines the recombination plot, gene bar and association scatter graph  
#' @param recombination.plot recombination plot 
#' @param gene.plot gene bar
#' @param marker.plot association scatter plot
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param ngenes number of genes in the genomic region
#' @param r2_legend add r2 legend
#' @import ggplot2 grid gridExtra gtable
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_assoc_combined3 <- function(recombination.plot, gene.plot, marker.plot, title=NULL, subtitle=NULL, ngenes, r2_legend=FALSE){
  #legend <- g_legend(marker.plot); marker.plot <- marker.plot + theme(legend.position="none")
  g1 <- ggplot_gtable(ggplot_build(recombination.plot))
  g2 <- ggplot_gtable(ggplot_build(marker.plot))
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g1$layout$name == "axis-l")
  ga <- g1$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  g$grobs[[3]] <- g2$grobs[[3]]
  g$grobs[[13]] <- g2$grobs[[13]]
  g <- gtable_add_cols(g, g1$widths[g1$layout[ia, ]$l], length(g$widths) - 1)
  if(ngenes<=10){
    g <- gtable_add_grob(g, list(textGrob("Recombination Rate (cM/Mb)", rot = -90, gp = gpar(col="black", fontsize=16))), pp$t, length(g$widths) - 1, pp$b)
  }else{
    g <- gtable_add_grob(g, list(textGrob("Recombination Rate", rot = -90, gp = gpar(col="black", fontsize=16))), pp$t, length(g$widths) - 1, pp$b)
  }
  g3 <- ggplot_gtable(ggplot_build(gene.plot))
  g3 <- gtable_add_grob(g3, g3$grobs[[which(g3$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  ia <- which(g3$layout$name == "axis-l")
  ga <- g3$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, ax, pp$t, length(g3$widths) - 1, pp$b)
  g3 <- gtable_add_cols(g3, g3$widths[g3$layout[ia, ]$l], length(g3$widths) - 1)
  g3 <- gtable_add_grob(g3, list(textGrob("", rot = -90, gp = gpar(fontsize=16, col = gray(.88)))), pp$t, length(g3$widths) - 1, pp$b)
  g <- gtable:::rbind_gtable(g, g3, "last")
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels[1]] <- unit(3,"null") 
  if(ngenes<=5){g$heights[panels[2]] <- unit(0.5,"null")}
  if(ngenes>5 & ngenes<=10){g$heights[panels[2]] <- unit(1.2,"null")}
  if(ngenes>10){g$heights[panels[2]] <- unit(2,"null")}
  if(!is.null(subtitle)){
    gt1 <- textGrob(subtitle,gp=gpar(fontsize=16))
    g <- gtable_add_rows(g, heights = grobHeight(gt1)*2.5, pos = 0)
    g <- gtable_add_grob(g, gt1, 1, 1, 1, ncol(g))
  }
  if(!is.null(title)){
    gt <- textGrob(title,gp=gpar(fontsize=20, fontface="italic"))
    g <- gtable_add_rows(g, heights = grobHeight(gt)*1.5, pos = 0)
    g <- gtable_add_grob(g, gt, 1, 1, 1, ncol(g))
  }
  g <- gtable_add_padding(g, unit(0.3, "cm"))
  if(r2_legend==TRUE){
    lheight <- sum(legend$height)*1.5
    g <- grid.arrange(g, legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
  }
  return(g)
}


##########################################################
##### Assoc plot #####
##########################################################

#' plot_assoc3
#'
#' plot_assoc3 plots a scatter graph of associations (e.g. log10 p-values) 
#' @param data data.frame with markername (marker), chromosome (chr), position (pos) and association statistics (stats)
#' @param corr correlation matrix between markers
#' @param corr.top correlation statistics between the top marker and the rest of the markers
#' @param x.min start of region
#' @param x.max end of region
#' @param top.marker the top associated marker, i.e. the marker with the largest -log10p or probability
#' @param ylab the y-axis label
#' @param type the type of the plot either log10p or probabilities
#' @import ggplot2
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
plot_assoc3 <- function(data, corr=NULL, corr.top=NULL, x.min, x.max, top.marker=NULL, ylab, type="log10p"){
  if(is.null(corr) & is.null(corr.top)) stop("no correlation statistics were input")
  if(is.null(corr) & !is.null(corr.top) & is.null(top.marker)) stop("top.marker must be defined if corr.top is provided")
  miss <- is.na(data$stats)
  if(!is.null(corr)){corr <- corr[!miss, !miss]}
  if(!is.null(corr.top)){corr.top <- corr.top[!miss]}  
  data <- data[!miss,]
  if(length(top.marker)!=0){
    top_marker <- which(top.marker==data$marker)
    if(length(top_marker)>1){top_marker <- sample(top_marker, 1); if(is.null(corr) & !is.null(corr.top)) warning("top.marker maps to multiple markers")}
    if(length(top_marker)==0){top_marker <- max.col(t(data$stats))} 
    lead_marker <- data[top_marker,]  
    ov_lead_marker <- data[max.col(t(data$stats)),]
    if((lead_marker$stats/ov_lead_marker$stats)>0.975){geomtext <- T}else{geomtext <- F}
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
  }else{ 
    top_marker <- max.col(t(data$stats))
    lead_marker <- data[top_marker,]  
    lead_marker$label_pos <- lead_marker$pos
    if((x.max-lead_marker$pos)<10000){lead_marker$label_pos <- lead_marker$pos - 0.025*(x.max-x.min)}
    if((lead_marker$pos-x.min)<10000){lead_marker$label_pos <- lead_marker$pos + 0.025*(x.max-x.min)}
    geomtext <- T
  }
  if(!is.null(corr)){r2 <- corr[,top_marker]^2}else{r2 <- corr.top^2}
  data$r2 <- "miss"
  data$r2[r2<0.2 & !is.na(r2)] <- "0.0-0.2"
  data$r2[r2>=0.2 & r2<=0.4 & !is.na(r2)] <- "0.2-0.4"
  data$r2[r2>=0.4 & r2<=0.6 & !is.na(r2)] <- "0.4-0.6"
  data$r2[r2>=0.6 & r2<=0.8 & !is.na(r2)] <- "0.6-0.8"
  data$r2[r2>=0.8 & r2<=1 & !is.na(r2)] <- "0.8-1.0" 
  data$r2 <- factor(data$r2, levels=c("miss", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
  ylim <- max(data$stats)
  marker.plot <- ggplot(aes(x=pos, y=stats), data=data) + geom_point(aes(fill=r2), pch=21, size=3.5) + scale_fill_manual(values=c("#DCDCDC", "#66FFFF", "#66FF66", "#FFCC00", "#FF9933", "#CC3300", "#FF0000"), drop=FALSE) + geom_point(data=lead_marker, aes(pos,stats), pch=22, colour="black", fill="black", size=4.2)  + theme_bw() + ylab(ylab) + xlab(NULL) + scale_y_continuous(limits=c(0,ylim)) + theme(axis.title.y=element_text(vjust=2.25, size=16), axis.text=element_text(size=14)) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_x_continuous(limits=c(x.min,x.max), breaks=NULL) + theme(axis.title=element_text(size=10)) + theme(legend.text=element_text(size=11), legend.title=element_text(size=12), legend.background = element_rect(colour = "black")) + theme(panel.background=element_rect(fill=NA)) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1)) 
  if(geomtext){marker.plot <- marker.plot + geom_text(data=lead_marker, aes(x=label_pos,y=stats,label=marker), vjust=-0.85, hjust=0.5, size=3.3)}else{if(lead_marker$stats[1]/ylim>=0.3){marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(-0.05*ylim), size=3.3, alpha=1)}else{marker.plot <- marker.plot + geom_label(data=lead_marker, aes(x=label_pos,y=stats,label=marker), label.r=unit(0, "lines"), nudge_y=(0.05*ylim), size=3.3, alpha=1)}} 
  if(type=="prob"){suppressMessages(marker.plot <- marker.plot + scale_y_continuous(limits=c(0,ylim), breaks=pretty(data$stats, n = 10)))}
  return(marker.plot)
}

