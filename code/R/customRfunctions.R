library(dendextend)

makePrettyMethCluster <- function(OBJ, metadata, my.cols.trt, my.cols.fam){
  ## Associate a color with a treatment
  metadata = join(metadata,
                  data.frame(trtG1G2 = unique(metadata$trtG1G2[order(metadata$trtG1G2)]),
                             coltrt = my.cols.trt))
  
  ## Associate a color with a family
  metadata = join(metadata,
                  data.frame(Family = unique(metadata$Family[order(metadata$Family)]),
                             colfam = my.cols.fam))
  ## Make dendrogram
  mydendro <- clusterSamples(OBJ, dist="correlation", method="ward", plot=FALSE)
  dend = as.dendrogram(mydendro)

  # Use color
  labels_colors(dend) <- metadata$coltrt[order.dendrogram(dend)]

  ## Hard coded methods here:
  dist.method="correlation"; hclust.method="ward.D"

  plot(dend, main = paste(OBJ@context, "methylation clustering"),
       sub = paste("Distance method: \"", dist.method,
                   "\"; Clustering method: \"", hclust.method,"\"",sep=""),
       xlab = "Samples", ylab = "Height")

  ## Ordered legend:
  correspDF=data.frame(name=unique(metadata$trtG1G2[order(metadata$trtG1G2)]),
                       color=my.cols.trt)
  legend("topright", title = "Treatment",
         legend = correspDF$name,
         col = correspDF$color, pch = 20,
         bty = "n",  pt.cex = 3, cex = 0.8 ,
         text.col = "black", horiz = FALSE, inset = c(0, 0.1))


  ## Add color bars and legend by family
  colored_bars(metadata$colfam, dend, # no need to order, it orders itself
               rowLabels = "Family")
  ## Legend
  legend("topleft", legend = unique(metadata$Family[mydendro$order]),
         fill = unique(metadata$colfam[mydendro$order]), bty = "n", title = "Family")
}
