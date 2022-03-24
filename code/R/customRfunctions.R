library(dendextend)
# for Manhattan plot:
library(dplyr)
library(tidyverse)

makePrettyMethCluster <- function(OBJ, metadata, my.cols.trt, my.cols.fam){
  ## Reorder metadata by sample ID, as OBJ methylkit!
  metadata = metadata[order(as.numeric(gsub("S", "", metadata$SampleID))),]
  
  ## Check
  if (!is.na(table(OBJ@sample.ids == metadata$SampleID)["FALSE"])){
    stop("check the samples order or similarity before both methylkit and metadata objects!")
  }
  
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
  # test
  print(metadata[metadata$SampleID %in% c("S37", "S41", "S46"), c("SampleID", "Family", "trtG1G2", "coltrt", "colfam")])
  # Use color
  labels_colors(dend) <- metadata$coltrt[order.dendrogram(dend)]

  ## Hard coded methods here:
  dist.method="correlation"; hclust.method="ward.D"

  plot(dend, main = paste(OBJ@context, "methylation clustering"),
       sub = paste("Distance method: \"", dist.method,
                   "\"; Clustering method: \"", hclust.method,"\"",sep=""),
       xlab = "Samples", ylab = "Height")

  ## Ordered legend:
  # correspDF=data.frame(name=unique(metadata$trtG1G2[order(metadata$trtG1G2)]),
  # correspDF=data.frame(name=unique(metadata$trtG1G2[order.dendrogram(dend)]),
  correspDF=data.frame(name=levels(metadata$trtG1G2),
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

makeManhattanPlots <- function(DMSfile, annotFile, GYgynogff, mycols=c("grey50","grey50","darkred","darkred"), 
                               mytitle = "Manhattan plot of DMS"){
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
  
  # table(DMSfile$chr)## check that chrXIX and chrUN are well removed!!
  
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
    scale_x_continuous(breaks=genome2$gmid,labels=genome2$chrom %>% str_remove(.,"Gy_chr"),
                       position = "top",expand = c(0,0))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line=element_blank(),
          axis.title = element_blank(),
          strip.placement = "outside")+
    ggtitle(mytitle)
}


###########################
## Adonis/NMDS functions ##
###########################
makePercentMetMat <- function(dataset){
  # creates a matrix containing percent methylation values
  perc.meth=percMethylation(dataset)
  # KOSTAS MBE: "Methylated sites and regions with low variation
  # and a standard deviation below 0.3, that is, noninformative
  # sites across individuals, were excluded from the cluster analyses"
  SD=apply(perc.meth,1, sd, na.rm = TRUE)
  perc.meth <- perc.meth[-which(SD<0.3),]
  x=t(perc.meth)
  return(x)
}

makeDatadistFUN <- function(dataset){
  x=makePercentMetMat(dataset)
  # creates a distance matrix. Method: Bray-Curtis, package vegan
  data.dist = as.matrix((vegdist(x, "bray", upper = FALSE))) 
}

myadonisFUN <- function(dataset, metadata){
  # make distance matrix with B-C distances
  data.dist = makeDatadistFUN(dataset)
  
  # We use a PERMANOVA to test the hypothesis that paternal treatment, 
  # offspring treatment, sex and their interactions significantly influencing global methylation
  perm <- how(nperm = 1000) # 1000 permutations
  setBlocks(perm) <- with(metadata, Family) # define the permutation structure considering family
  
  ## Full model
  print(adonis2(data.dist ~ PAT * outcome * Sex, data = metadata, permutations = perm))
  ## remove the non significant interactions
  print(adonis2(data.dist ~ PAT + outcome + Sex, data = metadata, permutations = perm))
  ## with the offspring treatments only (result of PAT + outcome)
  print(adonis2(data.dist ~ trtG1G2 + Sex, data = metadata, permutations = perm))
}

myGOF.NMDS.FUN <- function(dataset){
  # make distance matrix with B-C distances
  data.dist = makeDatadistFUN(dataset)
  # find the best number of dimensions (goeveg lib)
  ## Clarke 1993 suggests the following guidelines for acceptable stress values: <0.05 = excellent, <0.10
  # = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value
  # limit. Solutions with higher stress values should be interpreted with caution and those with stress
  # above 0.30 are highly suspect
  dimcheckMDS(
    data.dist,
    distance = "bray",
    k = 7,
    trymax = 100,
    autotransform = TRUE
  )
  abline(h = 0.1, col = "darkgreen")
}

myNMDS <- function(dataset, metadata){
  
  ## make percent methylation matrix
  x=makePercentMetMat(dataset)
  
  #Create NMDS based on bray-curtis distances - metaMDS finds the
  # most stable NMDS solution by randomly starting from different points in your data
  set.seed(38)
  NMDS <- metaMDS(comm = x, distance = "bray", maxit=1000, k = 6)
  
  #check to see stress of NMDS
  mystressplot <- stressplot(NMDS) 
  
  #extract plotting coordinates
  MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]
  ## OR #extract NMDS scores (x and y coordinates)
  ## data.scores = as.data.frame(scores(NMDS))
  
  #create new data table (important for later hulls finding)
  # with plotting coordinates and variables to test (dim 1,2,3)
  NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                   ID = metadata$ID,
                                   PAT=as.factor(metadata$PAT), 
                                   outcome=as.factor(metadata$outcome), 
                                   Sex = as.factor(metadata$Sex))
  #### start sub fun 
  makeNMDSplots <- function(dim, myvar){
    if (dim == "1_2"){
      dima=1; dimb=2
    } else if (dim == "1_3"){
      dima=1; dimb=3
    } else if (dim == "2_3"){
      dima=2; dimb=3
    }
    
    if (myvar == "PAT"){
      mycols = c("black","yellow")
    } else if (myvar == "Sex"){
      mycols = c("pink","blue")
    } else if (myvar == "outcome"){
      mycols = c("grey","red")
    }
    
    # generating convex hulls splitted by myvar in my metadata:
    hulls <- NMDS_dt[, .SD[chull(get(paste0("MDS", dima)), get(paste0("MDS", dimb)))], by = get(myvar)]
    
    myNMDSplot <- ggplot(NMDS_dt, 
                         aes_string(x=paste0("MDS",dima), y=paste0("MDS",dimb))) +
      geom_polygon(data = hulls, aes_string(fill=myvar), alpha=0.3) +
      scale_color_manual(values = mycols)+
      scale_fill_manual(values = mycols)+
      geom_point(aes_string(fill=myvar, shape=myvar), size = 3, alpha = .6) +
      #  geom_label(aes(label=row.names(NMDS2)))+
      scale_shape_manual(values = c(21,22)) +
      theme(legend.title=element_blank(), legend.position = "top")
    
    return(myNMDSplot)
  }
  
  figure <-  ggarrange(makeNMDSplots(dim= "1_2", myvar = "PAT"),
                       makeNMDSplots(dim= "1_3", myvar = "PAT"),
                       makeNMDSplots(dim= "2_3", myvar = "PAT"),
                       makeNMDSplots(dim= "1_2", myvar = "Sex"),
                       makeNMDSplots(dim= "1_3", myvar = "Sex"),
                       makeNMDSplots(dim= "2_3", myvar = "Sex"),
                       makeNMDSplots(dim= "1_2", myvar = "outcome"),
                       makeNMDSplots(dim= "1_3", myvar = "outcome"),
                       makeNMDSplots(dim= "2_3", myvar = "outcome"),
                       ncol = 3, nrow = 3)
  
  return(list(mystressplot=mystressplot, NMDSplot = figure))
}

## Calculate beta values (methylation proportion per CpG site) for the 1001880 positions covered in half G1 and half G2
getPMdataset <- function(uniteCov, MD, gener){
  PM = methylKit::percMethylation(uniteCov)
  
  ## Each row is a CpG sites, let's give them a proper "pos" row name
  rownames(PM) <- paste(uniteCov$chr, 
                        uniteCov$start, 
                        uniteCov$end)
  
  ## Select only the positions corresponding in DMS in G1 comparison control/infected
  length(DMS_info_G1$DMS)# 5074 DMS
  PM <- PM[rownames(PM) %in% DMS_info_G1$DMS, ]
  nrow(PM) # all good 
  
  ## Melt
  PM <- melt(PM)
  
  ## Extract chromosome, position, and assign correct names
  PM$Chr <- sapply(strsplit(as.character(PM$Var1), " +"), `[`, 1)
  PM$Pos <- sapply(strsplit(as.character(PM$Var1), " +"), `[`, 2)
  names(PM) <- c("Var1",  "ID",  "BetaValue", "Chr", "Pos")
  PM$rankpos <- 1:nrow(PM)
  
  ## Add treatment, Sex and Family
  dfTrt = data.frame(ID = MD$SampleID, Treatment = MD$trtG1G2, Sex = MD$Sex, Family= MD$Family)
  PM = merge(PM, dfTrt)
  
  if (gener=="parents"){
    PM$G1_trt <- PM$Treatment
    PM$G2_trt <- NA
  } else if (gener=="offspring"){
    PM$G1_trt <- sapply(strsplit(as.character(PM$Treatment), "_"), `[`, 1)
    PM$G2_trt <- sapply(strsplit(as.character(PM$Treatment), "_"), `[`, 2)
    PM$G1_trt[PM$G1_trt %in% "E"] <- "exposed"
    PM$G1_trt[PM$G1_trt %in% "NE"] <- "control"
  }
  ## Add the value of the DM in the parental comparison:
  names(PM)[names(PM) %in% "Var1"] <- "CpGSite"
  PM <- merge(PM, data.frame(CpGSite = DMS_info_G1$DMS, meth.diff.parentals = DMS_info_G1$meth.diff))
  
  ## Remove NA
  PM <- PM[!is.na(PM$BetaValue),]
  
  ## Add direction methylation diff in parental comparison
  PM$hypohyper <- "hypo"
  PM$hypohyper[PM$meth.diff.parentals > 0] <- "hyper"
  PM$hypohyper <- as.factor(PM$hypohyper)
  PM$hypohyper <- factor(PM$hypohyper, levels = c("hypo", "hyper"))
  
  return(PM)
}
