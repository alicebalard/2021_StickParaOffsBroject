makePrettyMethCluster <- function(OBJ, metadata){
  ## Clustering of raw CpG (those appearing on all samples):
  getColorVect<-function(x, colorPal){
    dfDat=data.frame(dat=x, num=as.numeric(x))
    dfCol=data.frame(color=colorPal, num=as.numeric(as.factor(colorPal)))
    df= join(dfDat,dfCol)
    df$color
  }
  
  # define palettes
  ## For Family, number can be 5 or 4
  if (length(unique(metadata$Family)==5)){
    seakingPal=c("#f69462","#5a5a62","#ffdebd","#943118","#cdc5c5")
  } else {    
    seakingPal=c("#f69462","#5a5a62","#ffdebd","#943118")
  }
  
  # Create a vector giving a color for each family (N=5 or 4)
  Family <- as.factor(metadata$Family)
  col_Family <- getColorVect(Family, seakingPal)
  
  ### plot
  mydendro <- clusterSamples(OBJ, dist="correlation", method="ward", plot=FALSE)
  
  # plot
  treatment=OBJ@treatment
  sample.ids=OBJ@sample.ids
  
  # order:"Control"    "E_control"  "E_exposed"  "Exposed"    "NE_control" "NE_exposed"
  treatmentName=unique(metadata$trtG1G2[order(metadata$trtG1G2)])
  my.cols=c("black", "blue", "purple","red","yellow","orange")
  
  col.list=as.list(my.cols[treatment])
  names(col.list)=sample.ids
  
  colLab <- function(n,col.list){
    if(is.leaf(n))
    {a <- attributes(n)
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col.list[[a$label]], lab.cex=1,
                                            col=col.list[[a$label]], cex=1, pch=16 ))
    }
    n
  }
  
  dend = as.dendrogram(mydendro)
  dend_colored <- dendrapply(dend, colLab,col.list)
  ## Hard coded methods here:
  dist.method="correlation"; hclust.method="ward.D"
  
  plot(dend_colored, main = paste(OBJ@context, "methylation clustering"), 
       sub = paste("Distance method: \"", dist.method,
                   "\"; Clustering method: \"", hclust.method,"\"",sep=""), 
       xlab = "Samples", ylab = "Height")
  
  ## To order legend:
  correspDF=data.frame(name=unique(treatmentName[order(treatmentName)]),
                       color=my.cols)
  correspDF=correspDF[c(1,4,5,6,2,3),]
  
  legend("topright", title = "Treatment",
         legend = correspDF$name,
         col = correspDF$color, pch = 20,
         bty = "n",  pt.cex = 3, cex = 0.8 , 
         text.col = "black", horiz = FALSE, inset = c(0, 0.1))
  
  ## Add color bars and legend
  colored_bars(col_Family, mydendro,
               rowLabels = "Family")
  
  ## Legend
  legend("topleft", legend = unique(metadata$Family[mydendro$order]),
         fill = unique(col_Family[mydendro$order]), bty = "n", title = "Family")
}

