library(dendextend)
# for Manhattan plot:
library(dplyr)
library(tidyverse)

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

####################
## Manhattan plot ##
####################
# Adapted from Melanie Heckwolf
# create Manhattan plots over the genome from a DMS file

#the stickleback genome is annotated as chrI, chrII, chrIII ... to chrUn (unknown, includes all scaffolds)
#this function takes the roman rumbers and turns them into actual numbers:
deroman <- function(x){ x %>% str_remove(.,"Gy_chr") %>% 
    ifelse(. %in% c("Un","M"),., as.roman(.) %>% 
             as.integer() %>% as.character() %>% str_pad(.,2,"left",0))
}

get_pos <- function(CHROM,genome){
  tibble(CHROM = CHROM, 
         POS = sample(x=1:genome$length[genome$chrom == CHROM],size=1))
}

makeManhattanPlots <- function(DMSfile, annotFile, GYgynogff, mycols=c("grey50","grey50","darkred","darkred")){
  #GA_genome.fa.sizes.txt is a file with chromosome sizes and names
  genome <- GYgynogff %>%
    mutate(chrom_nr=chrom %>% deroman(),
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
    arrange(chrom_order) %>%
    mutate(gstart=lag(length,default=0) %>% cumsum(),
           gend=gstart+length, 
           type=LETTERS[2-(chrom_order%%2)],
           gmid=(gstart+gend)/2)
  
  #genome without M re-type:
  genome2=genome[genome$chrom_nr!="M",] %>%
    mutate(type=rep(c("A","B"),length(length)/2))
  
  region=as.factor(ifelse(annotFile$prom==1,"promoter",
                          ifelse(annotFile$exon==1,"exon",
                                 ifelse(annotFile$intron==1, "intron","intergenic"))))
  
  mydata = tibble(chrom=DMSfile$chr,
                  pos=DMSfile$start,
                  meth.diff=DMSfile$meth.diff,
                  qval=DMSfile$qvalue,
                  region=region)
  
  table(DMSfile$chr)## check that chrXIX and chrUN are well removed!!
  
  # join DMS and genomic position
  data = left_join(mydata, genome2) %>% 
    mutate(gpos=pos+gstart,significance= ifelse(abs(qval>0.0125) | abs(meth.diff)<15,"not significant","significant"))
  
  table(data$significance) # all signif
  
  #plot only significant DMS:
  ggplot()+
    geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
    geom_point(data=data[abs(data$meth.diff)>15 & data$significance=="significant",],
               aes(x=gpos,y=meth.diff,col=region,shape=region),fill="white", size = 2)+
    scale_color_manual(values = mycols)+
    scale_shape_manual(values=c(21,21,21,21))+
    scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
    scale_x_continuous(breaks=genome2$gmid,labels=genome2$chrom %>% str_remove(.,"chr"),
                       position = "top",expand = c(0,0))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line=element_blank(),
          axis.title = element_blank(),
          strip.placement = "outside")
}