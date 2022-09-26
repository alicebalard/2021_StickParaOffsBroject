## I. Functions used in R03.2

#########################
## Clustering function ##
#########################

makePrettyMethCluster <- function(OBJ, metadata, my.cols.trt, my.cols.fam, nbrk){
  ## Reorder metadata by sample ID, as OBJ methylkit!
  metadata = metadata[order(as.numeric(gsub("S", "", metadata$SampleID))),]
  
  ## Check
  if (!is.na(table(OBJ@sample.ids == metadata$SampleID)["FALSE"])){
    stop("check the samples order or similarity before both methylkit and metadata objects!")
  }
  
  ## To add color bars:
  # Generate color palette
  pal = qualpalr::qualpal(55, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))
  ## Family
  fam <- factor(metadata$Family)
  col_fam <- pal$hex[1:length(unique(metadata$Family))][fam]
  ## Treatment
  trt <- factor(metadata$outcome)
  col_trt <- c("grey", "red")[trt]
  ## Paternal treatment
  trtPAT <- factor(metadata$PAT)
  col_trtPAT <- c("grey", "red")[trtPAT]
  ## Brother pair of the father
  brotherPairID <- factor(metadata$brotherPairID)
  x = length(levels(brotherPairID))
  col_brotherPairID <- sample(pal$hex, x)[brotherPairID]
  ## Clutch ID
  # clutch <- factor(metadata$clutch.ID)
  # x = length(levels(clutch))
  # col_clutch <- sample(pal$hex, x)[clutch]
  
  ## Make dendrogram
  mydendro <- clusterSamples(OBJ, dist="correlation", method="ward", plot=FALSE)
  dend = as.dendrogram(mydendro)
  
  ## and plot
  dend %>% plot(main=paste(OBJ@context, "methylation clustering\n", 
                           "Distance method: correlation; Clustering method: ward.D"), 
                font.main = 1, cex.main = .8, ylab = "Height", nodePar = list(lab.cex = 0.6, pch = c(NA, NA))) 
  dend %>% rect.dendrogram(k=nbrk, border = 8, lty = 5, lwd = 2)
  colored_bars(cbind(col_trt, col_trtPAT, brotherPairID, col_fam), dend, y_shift = -0.1, #col_clutch
               rowLabels = c("G2 treatment", "G1 treatment", "G1 family",  "G0 family")) #"Clutch"
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

######################
## Adonis functions ##
######################
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

AdonisWithinG1trtFUN <- function(trtgp){
  # make distance matrix with B-C distances
  data.dist = makeDatadistFUN(reorganize(methylObj = uniteCovALL_G2_woSexAndUnknowChr,
                                         treatment = fullMetadata_OFFS$trtG1G2_NUM[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp],
                                         sample.ids = fullMetadata_OFFS$ID[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp]))
  perm <- how(nperm = 1000) # 1000 permutations
  # define the permutation structure considering brotherPairID
  setBlocks(perm) <- with(fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], brotherPairID)
  adonis2(data.dist ~ outcome + Sex + brotherPairID,
          data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2_NUM %in% trtgp,], permutations = perm)
}

####################
## NMDS functions ##
####################
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

myNMDSFUN <- function(dataset, metadata, myseed, byParentTrt=FALSE, trtgp=NA){
  
  print(paste0("my seed = ", myseed))
  
  if (byParentTrt==TRUE){
    dataset = reorganize(methylObj = dataset,
                         treatment = metadata$trtG1G2_NUM[metadata$trtG1G2_NUM %in% trtgp],
                         sample.ids = metadata$ID[metadata$trtG1G2_NUM %in% trtgp])
    metadata = metadata[metadata$trtG1G2_NUM %in% trtgp, ]
  }
  
  ## make percent methylation matrix
  x=makePercentMetMat(dataset)
  #Create NMDS based on bray-curtis distances - metaMDS finds the
  # most stable NMDS solution by randomly starting from different points in your data
  set.seed(myseed)
  NMDS <- metaMDS(comm = x, distance = "bray", maxit=1000, k = 6)
  #check to see stress of NMDS
  mystressplot <- stressplot(NMDS) 
  #extract plotting coordinates
  MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]
  ## OR #extract NMDS scores (x and y coordinates)
  ## data.scores = as.data.frame(scores(NMDS))
  
  #create new data table (important for later hulls finding)
  # with plotting coordinates and variables to test (dim 1,2,3)
  
  if (byParentTrt==FALSE){
    NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                     ID = metadata$ID,
                                     PAT=as.factor(metadata$PAT), 
                                     outcome=as.factor(metadata$outcome), 
                                     Sex = as.factor(metadata$Sex),
                                     brotherPairID = as.factor(metadata$brotherPairID))
  } else if (byParentTrt==TRUE){
    NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                     ID = metadata$ID,
                                     outcome=as.factor(metadata$outcome), 
                                     Sex = as.factor(metadata$Sex),
                                     brotherPairID = as.factor(metadata$brotherPairID))
  }
  
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
      mycols = c("black","yellow"); myshape = c(21,22)
    } else if (myvar == "Sex"){
      mycols = c("pink","blue"); myshape = c(21,22)
    } else if (myvar == "outcome"){
      mycols = c("grey","red"); myshape = c(21,22)
    } else if (myvar == "brotherPairID"){
      mycols = c(1:8); myshape = rep(21, 8)
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
      scale_shape_manual(values = myshape) +
      theme(legend.title=element_blank(), legend.position = "top")
    
    return(myNMDSplot)
  }
  
  if (byParentTrt==FALSE){
    figure <-  ggarrange(makeNMDSplots(dim= "1_2", myvar = "PAT"),
                         makeNMDSplots(dim= "1_3", myvar = "PAT"),
                         makeNMDSplots(dim= "2_3", myvar = "PAT"),
                         makeNMDSplots(dim= "1_2", myvar = "Sex"),
                         makeNMDSplots(dim= "1_3", myvar = "Sex"),
                         makeNMDSplots(dim= "2_3", myvar = "Sex"),
                         makeNMDSplots(dim= "1_2", myvar = "outcome"),
                         makeNMDSplots(dim= "1_3", myvar = "outcome"),
                         makeNMDSplots(dim= "2_3", myvar = "outcome"),
                         makeNMDSplots(dim= "1_2", myvar = "brotherPairID"),
                         makeNMDSplots(dim= "1_3", myvar = "brotherPairID"),
                         makeNMDSplots(dim= "2_3", myvar = "brotherPairID"),
                         ncol = 3, nrow = 4)
  } else if (byParentTrt==TRUE){
    figure <-  ggarrange(makeNMDSplots(dim= "1_2", myvar = "Sex"),
                         makeNMDSplots(dim= "1_3", myvar = "Sex"),
                         makeNMDSplots(dim= "2_3", myvar = "Sex"),
                         makeNMDSplots(dim= "1_2", myvar = "outcome"),
                         makeNMDSplots(dim= "1_3", myvar = "outcome"),
                         makeNMDSplots(dim= "2_3", myvar = "outcome"),
                         makeNMDSplots(dim= "1_2", myvar = "brotherPairID"),
                         makeNMDSplots(dim= "1_3", myvar = "brotherPairID"),
                         makeNMDSplots(dim= "2_3", myvar = "brotherPairID"),
                         ncol = 3, nrow = 3)
  }
  return(list(NMDS = NMDS, mystressplot=mystressplot, NMDSplot = figure))
}

#########
## PCA ## 
#########
myPCA <- function(x, incomplete){
  if (incomplete==TRUE){
    # estimate the number of components from incomplete data
    nb <- estim_ncpPCA(x, scale = T)
    # impute the table
    res.comp <- imputePCA(x, ncp = nb$ncp, scale = T)
    x = res.comp$completeObs
  }
  # 2. run PCA
  res.PCA = FactoMineR::PCA(x, scale.unit = T, graph = FALSE) # perform PCA
  metadata = fullMetadata_OFFS
  # check that the sample names are in the same order
  ifelse(table(rownames(res.PCA$ind$coord) == metadata$SampleID), "sample names are in the same order", "ERROR in PCA sample names order")
  # 3. extract axes 1, 2
  metadata$PCA1 = res.PCA$ind$coord[,1] # axis 1
  metadata$PCA2 = res.PCA$ind$coord[,2] # axis 2
  # 4. Correlation with parasite load/BCI
  mod = lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
  ## Model selection:
  modSel = lmer(formula = attr(attr(lmerTest::step(mod, reduce.random = F), "drop1"), "heading")[3],
                data=metadata, REML = F)
  print("The chosen model is:")
  print(formula(modSel))
  return(list(res.PCA=res.PCA, modSel = modSel, metadata = metadata))
}

getPCACpG <- function(DMSvec, effect){
  pos2keep = which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, uniteCov14_G2_woSexAndUnknowChrOVERLAP$start, sep = " ") %in%
                     DMSvec)
  uniteAtDMS = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, pos2keep)
  percAtDMS = percMethylation(uniteAtDMS)
  
  print(paste(nrow(percAtDMS), "DMS linked with", effect))
  
  # We use missMDA and FactoMineR for imputation of missing data and
  # performing of PCA: (see:
  #
  #                       -   <http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html>
  #                       -   <https://www.youtube.com/watch?v=OOM8_FH6_8o>)
  
  PCA_percAtDMS_imputed <- myPCA(x = t(percAtDMS), incomplete = TRUE)
  
  # The function dimdesc() can be used to identify the most correlated variables with a given principal component.
  mydimdesc = dimdesc(PCA_percAtDMS_imputed$res.PCA, axes = c(1,2), proba = 0.05)
  
  print(paste(nrow(mydimdesc$Dim.1$quanti), "CpG sites most correlated (p < 0.05) with the first principal component"))
  print(paste(nrow(mydimdesc$Dim.2$quanti), "CpG sites most correlated (p < 0.05) with the second principal component"))
  
  # Extract the values for CpGs associated with the
  CpGPCA1 = methylKit::select(uniteAtDMS, as.numeric(gsub("V","",rownames(mydimdesc$Dim.1$quanti))))
  CpGPCA2 = methylKit::select(uniteAtDMS, as.numeric(gsub("V","",rownames(mydimdesc$Dim.2$quanti))))
  
  return(list(PCA_percAtDMS_imputed=PCA_percAtDMS_imputed, CpGPCA1=CpGPCA1, CpGPCA2=CpGPCA2))
}

##########################
## Calculate beta values (methylation proportion per CpG site) for the 1001880 positions covered in half G1 and half G2
getPMdataset <- function(uniteCov, MD, gener){
  PM = methylKit::percMethylation(uniteCov)
  
  ## Each row is a CpG sites, let's give them a proper "pos" row name
  rownames(PM) <- paste(uniteCov$chr, uniteCov$start, uniteCov$end)
  
  ## Select only the positions corresponding in DMS in G1 comparison control/infected
  length(DMS_info_G1$DMS)
  PM <- PM[rownames(PM) %in% DMS_info_G1$DMS, ]
  nrow(PM) # all good 
  
  ## Melt
  PM <- melt(PM)
  
  ## Extract chromosome, position, and assign correct names
  PM$Chr <- sapply(strsplit(as.character(PM$Var1), " +"), `[`, 1)
  PM$Pos <- sapply(strsplit(as.character(PM$Var1), " +"), `[`, 2)
  names(PM) <- c("Var1",  "ID",  "BetaValue", "Chr", "Pos")
  PM$rankpos <- 1:nrow(PM)
  
  ## Add treatment, Sex, brotherPairID and clutchID
  dfTrt = data.frame(ID = MD$SampleID, Treatment = MD$trtG1G2, Sex = MD$Sex, brotherPairID= MD$brotherPairID, clutch.ID = MD$clutch.ID)
  PM = merge(PM, dfTrt)
  
  if (gener=="parents"){
    PM$G1_trt <- PM$Treatment
    PM$G2_trt <- NA
  } else if (gener=="offspring"){
    PM$G1_trt <- sapply(strsplit(as.character(PM$Treatment), "_"), `[`, 1)
    PM$G2_trt <- sapply(strsplit(as.character(PM$Treatment), "_"), `[`, 2)
    PM$G2_trt[PM$G2_trt %in% "control"] <- "Control"
    PM$G2_trt[PM$G2_trt %in% "exposed"] <- "Exposed"
    PM$G1_trt[PM$G1_trt %in% "E"] <- "Exposed"
    PM$G1_trt[PM$G1_trt %in% "NE"] <- "Control"
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

## Additional function for complexUpset to color by degrees
query_by_degree = function(data, groups, params_by_degree, ...) {
  intersections = unique(ComplexUpset::upset_data(data, groups)$plot_intersections_subset)
  lapply(
    intersections,
    FUN=function(x) {
      members = strsplit(x, '-', fixed=TRUE)[[1]]
      if (!(length(members) %in% names(params_by_degree))) {
        stop(
          paste('Missing specification of params for degree', length(members))
        )
      }
      args = c(
        list(intersect=members, ...),
        params_by_degree[[length(members)]]
      )
      do.call(ComplexUpset::upset_query, args)
    }
  )
}

########################################
## Differential methylation functions ##
########################################
getDiffMeth <- function(myuniteCov, myMetadata, mccores=10, mydif = 15){
  if (length(table(myMetadata$Sex)) == 1 & length(table(myMetadata$brotherPairID)) == 1){ # 1 sex, 1 BP -> no covariate
    myDiffMeth=calculateDiffMeth(myuniteCov, mc.cores = mccores)#10 on Apocrita
  } else { # if more than 1 sex or 1 BP, we add a covariate
    if (length(table(myMetadata$Sex)) == 1 & length(table(myMetadata$brotherPairID)) > 1){
      cov = data.frame(brotherPairID = myMetadata$brotherPairID)
    } else if (length(table(myMetadata$Sex)) == 2 & length(table(myMetadata$brotherPairID)) > 1){
      cov = data.frame(brotherPairID = myMetadata$brotherPairID, Sex = myMetadata$Sex)
    } else if (length(table(myMetadata$Sex)) == 2){ # this is for within brother pairs
      cov = data.frame(Sex = myMetadata$Sex)
    } 
    myDiffMeth=calculateDiffMeth(myuniteCov, covariates = cov, mc.cores = mccores)#10 on Apocrita
  }
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
  ## NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
  myDMS_15pc = getMethylDiff(myDiffMeth, difference=mydif, qvalue=0.01)
  return(myDMS_15pc)
}

getDiffMethSimple <- function(myuniteCov, myMetadata){
  myDiffMeth=calculateDiffMeth(myuniteCov, mc.cores = 3)#10 on Apocrita
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%.
  ## NB: arg type="hyper" or type="hypo" gives hyper-methylated or hypo-methylated regions/bases.
  myDMS_15pc = getMethylDiff(myDiffMeth, difference=15, qvalue=0.01)
  return(myDMS_15pc)
}

## Calculate DMS in all brother pairs
## 1. CC-TC = CONTROL fish (parent CvsT)
## 2. CT-TT = TREATMENT fish (parent CvsT)
## 3. CC-CT = fish from CONTROL parents (G2 CvsT)
## 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

## Per brother pair:
getDMperBP <- function(BP){
  ## Unite object for one Brother Pair:
  metadataBP_CC_TC = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_control", "E_control"), ]
  metadataBP_CT_TT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"), ]
  metadataBP_CC_CT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("NE_control", "NE_exposed"), ]
  metadataBP_TC_TT = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                         fullMetadata_OFFS$trtG1G2 %in% c("E_control", "E_exposed"), ]
  
  ## Make 4 separate uniteCov:
  myuniteCovBP_CC_TC = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CC_TC$trtG1G2_NUM, sample.ids = metadataBP_CC_TC$ID)
  myuniteCovBP_CT_TT = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CT_TT$trtG1G2_NUM, sample.ids = metadataBP_CT_TT$ID)
  myuniteCovBP_CC_CT = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_CC_CT$trtG1G2_NUM, sample.ids = metadataBP_CC_CT$ID)
  myuniteCovBP_TC_TT = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                  treatment = metadataBP_TC_TT$trtG1G2_NUM, sample.ids = metadataBP_TC_TT$ID)
  
  ## remove bases where NO fish in this BP has a coverage
  myuniteCovBP_CC_TC = methylKit::select(myuniteCovBP_CC_TC, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_TC)))))
  myuniteCovBP_CT_TT = methylKit::select(myuniteCovBP_CT_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CT_TT)))))
  myuniteCovBP_CC_CT = methylKit::select(myuniteCovBP_CC_CT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_CT)))))
  myuniteCovBP_TC_TT = methylKit::select(myuniteCovBP_TC_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_TC_TT)))))
  
  ## Calculate differential methylation:
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate
  DMS_15pc_BP_CC_TC = getDiffMeth(myuniteCov = myuniteCovBP_CC_TC, myMetadata = metadataBP_CC_TC, mccores = 10, mydif = 15)
  DMS_15pc_BP_CT_TT = getDiffMeth(myuniteCov = myuniteCovBP_CT_TT, myMetadata = metadataBP_CT_TT, mccores = 10, mydif = 15)
  DMS_15pc_BP_CC_CT = getDiffMeth(myuniteCov = myuniteCovBP_CC_CT, myMetadata = metadataBP_CC_CT, mccores = 10, mydif = 15)
  DMS_15pc_BP_TC_TT = getDiffMeth(myuniteCov = myuniteCovBP_TC_TT, myMetadata = metadataBP_TC_TT, mccores = 10, mydif = 15)
  
  ## tile for Differentially methylated REGIONS
  tilesBP_CC_TC = tileMethylCounts(myuniteCovBP_CC_TC,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_TC = getDiffMeth(tilesBP_CC_TC, metadataBP_CC_TC)
  tilesBP_CT_TT = tileMethylCounts(myuniteCovBP_CT_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CT_TT = getDiffMeth(tilesBP_CT_TT, metadataBP_CT_TT)
  tilesBP_CC_CT = tileMethylCounts(myuniteCovBP_CC_CT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_CT = getDiffMeth(tilesBP_CC_CT, metadataBP_CC_CT)
  tilesBP_TC_TT = tileMethylCounts(myuniteCovBP_TC_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_TC_TT = getDiffMeth(tilesBP_TC_TT, metadataBP_TC_TT)
  
  return(list(DMSlist = list(DMS_15pc_BP_CC_TC = DMS_15pc_BP_CC_TC, DMS_15pc_BP_CT_TT = DMS_15pc_BP_CT_TT, DMS_15pc_BP_CC_CT = DMS_15pc_BP_CC_CT, DMS_15pc_BP_TC_TT = DMS_15pc_BP_TC_TT),
              DMRlist = list(DMR_15pc_BP_CC_TC = DMR_15pc_BP_CC_TC, DMR_15pc_BP_CT_TT = DMR_15pc_BP_CT_TT, DMR_15pc_BP_CC_CT = DMR_15pc_BP_CC_CT, DMR_15pc_BP_TC_TT = DMR_15pc_BP_TC_TT)))
}

## And for positions covered in ALL FISH
getDMperBP2 <- function(BP){
  ## Unite object for one Brother Pair:
  metadataBP_CC_TC = fullMetadata[fullMetadata$brotherPairID %in% BP &
                                    fullMetadata$trtG1G2 %in% c("NE_control", "E_control"), ]
  metadataBP_CT_TT = fullMetadata[fullMetadata$brotherPairID %in% BP &
                                    fullMetadata$trtG1G2 %in% c("NE_exposed", "E_exposed"), ]
  metadataBP_CC_CT = fullMetadata[fullMetadata$brotherPairID %in% BP &
                                    fullMetadata$trtG1G2 %in% c("NE_control", "NE_exposed"), ]
  metadataBP_TC_TT = fullMetadata[fullMetadata$brotherPairID %in% BP &
                                    fullMetadata$trtG1G2 %in% c("E_control", "E_exposed"), ]
  
  ## Make 4 separate uniteCov:
  myuniteCovBP_CC_TC = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                                  treatment = metadataBP_CC_TC$trtG1G2_NUM, sample.ids = metadataBP_CC_TC$ID)
  myuniteCovBP_CT_TT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                                  treatment = metadataBP_CT_TT$trtG1G2_NUM, sample.ids = metadataBP_CT_TT$ID)
  myuniteCovBP_CC_CT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                                  treatment = metadataBP_CC_CT$trtG1G2_NUM, sample.ids = metadataBP_CC_CT$ID)
  myuniteCovBP_TC_TT = reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                                  treatment = metadataBP_TC_TT$trtG1G2_NUM, sample.ids = metadataBP_TC_TT$ID)
  
  ## remove bases where NO fish in this BP has a coverage
  myuniteCovBP_CC_TC = methylKit::select(myuniteCovBP_CC_TC, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_TC)))))
  myuniteCovBP_CT_TT = methylKit::select(myuniteCovBP_CT_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CT_TT)))))
  myuniteCovBP_CC_CT = methylKit::select(myuniteCovBP_CC_CT, which(!is.na(rowSums(percMethylation(myuniteCovBP_CC_CT)))))
  myuniteCovBP_TC_TT = methylKit::select(myuniteCovBP_TC_TT, which(!is.na(rowSums(percMethylation(myuniteCovBP_TC_TT)))))
  
  ## Calculate differential methylation:
  ## We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate
  DMS_15pc_BP_CC_TC = getDiffMeth(myuniteCov = myuniteCovBP_CC_TC, myMetadata = metadataBP_CC_TC, mccores = 10, mydif = 15)
  DMS_15pc_BP_CT_TT = getDiffMeth(myuniteCov = myuniteCovBP_CT_TT, myMetadata = metadataBP_CT_TT, mccores = 10, mydif = 15)
  DMS_15pc_BP_CC_CT = getDiffMeth(myuniteCov = myuniteCovBP_CC_CT, myMetadata = metadataBP_CC_CT, mccores = 10, mydif = 15)
  DMS_15pc_BP_TC_TT = getDiffMeth(myuniteCov = myuniteCovBP_TC_TT, myMetadata = metadataBP_TC_TT, mccores = 10, mydif = 15)
  
  ## tile for Differentially methylated REGIONS
  tilesBP_CC_TC = tileMethylCounts(myuniteCovBP_CC_TC,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_TC = getDiffMeth(tilesBP_CC_TC, metadataBP_CC_TC)
  tilesBP_CT_TT = tileMethylCounts(myuniteCovBP_CT_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CT_TT = getDiffMeth(tilesBP_CT_TT, metadataBP_CT_TT)
  tilesBP_CC_CT = tileMethylCounts(myuniteCovBP_CC_CT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_CC_CT = getDiffMeth(tilesBP_CC_CT, metadataBP_CC_CT)
  tilesBP_TC_TT = tileMethylCounts(myuniteCovBP_TC_TT,win.size=100,step.size=100,cov.bases = 10)
  DMR_15pc_BP_TC_TT = getDiffMeth(tilesBP_TC_TT, metadataBP_TC_TT)
  
  return(list(DMSlist = list(DMS_15pc_BP_CC_TC = DMS_15pc_BP_CC_TC, DMS_15pc_BP_CT_TT = DMS_15pc_BP_CT_TT, DMS_15pc_BP_CC_CT = DMS_15pc_BP_CC_CT, DMS_15pc_BP_TC_TT = DMS_15pc_BP_TC_TT),
              DMRlist = list(DMR_15pc_BP_CC_TC = DMR_15pc_BP_CC_TC, DMR_15pc_BP_CT_TT = DMR_15pc_BP_CT_TT, DMR_15pc_BP_CC_CT = DMR_15pc_BP_CC_CT, DMR_15pc_BP_TC_TT = DMR_15pc_BP_TC_TT)))
}

## Function to get DMS info
myDMSinfo <- function(DMSobject, fromUniteCov){
  DMS = paste(DMSobject$chr, DMSobject$start, DMSobject$end)
  meth.diff = DMSobject$meth.diff
  direction = ifelse(DMSobject$meth.diff > 0, "hyper", "hypo")
  percentDMS = length(DMS)/nrow(fromUniteCov)*100
  return(list(DMS = DMS, meth.diff = meth.diff, direction = direction, percentDMS = percentDMS))
}

# order positions by chromosomes & position
reorderByChrom <- function(x){
  df = data.frame(fullpos=names(x), beta=x, row.names = NULL)
  df$chr = sapply(strsplit(df$fullpos,"\\."), `[`, 1)
  df$pos = sapply(strsplit(df$fullpos,"\\."), `[`, 2)
  df = df %>%
    mutate(chrom_nr=chr %>% deroman(), # deroman is custom, defined in customRfunctions.R
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
    arrange(chrom_order) 
  orderedVec = df$beta
  names(orderedVec) = df$fullpos
  return(orderedVec)
}

# calculate average methylation per treatment group at each position
calcAveMeth <- function(perc_uniteObj){
  rawmetadata = fullMetadata[match(colnames(perc_uniteObj), fullMetadata$SampleID), ]
  
  for (i in 1:length(levels(rawmetadata$trtG1G2))){
    whichCols = which(colnames(perc_uniteObj) %in% rawmetadata$SampleID[
      rawmetadata$trtG1G2 %in% levels(rawmetadata$trtG1G2)[i]])
    perc_uniteObj = perc_uniteObj %>% data.frame() %>%
      dplyr::mutate(X = rowMeans(dplyr::select(., whichCols), na.rm = T))
    names(perc_uniteObj)[names(perc_uniteObj) %in% "X"] = paste0("ave", levels(rawmetadata$trtG1G2)[i])
  }
  perc_uniteObj = perc_uniteObj[grep("ave", names(perc_uniteObj))]
}

## Manhattan plots function:
# GYgynogff a data frame with a "chrom" and a "length" columns 
# (NB: here "genome4Manhattan" is specific to my stickleback file)
plotManhattanGenesDMS <- function(annotFile, GYgynogff){
  annotFile=annotFile %>% 
    dplyr::select(c("start.gene", "end.gene", "GeneSymbol", "feature.name", "Note", "chrom",
                    "nDMSperGenekb", "ENTREZID", "description", "summary"))%>% unique
  
  ## Prepare genome for Manhattan plots:
  genome4Manhattan = GYgynogff %>%
    #genome without chrXIX and unknown re-type:
    filter(chrom!="Gy_chrXIX" & chrom!= "Gy_chrUn")%>%
    mutate(chrom_nr=chrom %>% deroman(), 
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% arrange(chrom_order) %>%
    mutate(gstart=lag(length,default=0) %>% cumsum(), 
           gend=gstart+length, 
           typeBG=LETTERS[2-(chrom_order%%2)],   
           gmid=(gstart+gend)/2)
  
  # Prepare data and change gene position to start at the good chromosome
  data4Manhattan = dplyr::left_join(annotFile, genome4Manhattan) %>% 
    dplyr::mutate(posInPlot=((end.gene+start.gene)/2)+gstart)
  
  # Short name of the gene if we want to plot these as labels:
  data4Manhattan$Note = unlist(data4Manhattan$Note)
  data4Manhattan$Note = str_extract(data4Manhattan$Note, "(?<=Similar to )(\\w+)")
  
  # Manhattan plot
  plot = ggplot()+
    # add grey background every second chromosome
    geom_rect(data=genome4Manhattan,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=typeBG), alpha=.2)+
    scale_x_continuous(breaks=genome4Manhattan$gmid,labels=genome4Manhattan$chrom %>% str_remove(.,"Gy_chr"),
                       position = "top",expand = c(0,0))+
    scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none") +
    # geom_hline(yintercept = 1)+ # if want to add line break
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ # add frame
    # scale_y_continuous(breaks = seq(0, 20), expand = expansion(mult = 0.5)) + # increase size under plot for labels
    ylab("Number of differentially methylated CpG per gene kb")+ 
    geom_point(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGenekb)) +
    geom_label_repel(data = data4Manhattan[data4Manhattan$nDMSperGenekb > 1,],
                     aes(x=posInPlot, y = nDMSperGenekb, label = Note), max.overlaps = Inf)
  return(plot)
}

# GYgynogff a data frame with a "chrom" and a "length" columns (NB: here "genome4Manhattan" is specific to my stickleback file)
plotManhattanGenesDMS4BP <- function(annotFile, i = 0, GYgynogff, myxlab = NULL, isBPinfo=TRUE){
  ## Prepare genome for Manhattan plots:
  genome4Manhattan = GYgynogff %>%
    #genome without chrXIX and unknown re-type:
    filter(chrom!="Gy_chrXIX" & chrom!= "Gy_chrUn")%>%
    mutate(chrom_nr=chrom %>% deroman(), 
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% arrange(chrom_order) %>%
    mutate(gstart=lag(length,default=0) %>% cumsum(), 
           gend=gstart+length, 
           typeBG=LETTERS[2-(chrom_order%%2)],   
           gmid=(gstart+gend)/2)
  
  # Prepare data and change gene position to start at the good chromosome
  data4Manhattan = dplyr::left_join(annotFile, genome4Manhattan) %>% dplyr::mutate(posInPlot=start+gstart)
  
  # Short name of the gene if we want to plot these as labels:
  data4Manhattan$Note = unlist(data4Manhattan$Note)
  data4Manhattan$Note = str_extract(data4Manhattan$Note, "(?<=Similar to )(\\w+)")
  
  # Manhattan plot
  plot = ggplot()+
    # add grey background every second chromosome
    geom_rect(data=genome4Manhattan,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=typeBG), alpha=.2)+
    scale_x_continuous(breaks=genome4Manhattan$gmid,labels=genome4Manhattan$chrom %>% str_remove(.,"Gy_chr"),
                       position = "top",expand = c(0,0))+
    scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none") +
    # geom_hline(yintercept = 1)+ # if want to add line break
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ # add frame
    scale_y_continuous(breaks = seq(0, 20), expand = expansion(mult = 0.5)) + # increase size under plot for labels
    ylab("Number of differentially methylated CpG per gene kb")
  
  # add points
  if (isBPinfo==TRUE){
    plot = plot + 
      geom_point(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGenekb, col=as.factor(nbrBP)), size = 2) +
      scale_color_manual(values = c('grey', 'red', 'purple', 'blue', 'green'),
                         name = "Genes found differentially methylated in N brother pairs:") +
      xlab(paste0("Genes with DMS present in at least 4 brother pairs\nComparison: ", vecCompa[i]))
  } else {
    plot = plot + 
      geom_point(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGenekb), size = 2) +
      geom_label_repel(data = data4Manhattan, aes(x=posInPlot, y = nDMSperGenekb, label = Note), max.overlaps = Inf)+
      xlab(myxlab)
  }
  return(plot)
}

########
## GO ## 
########
makedfGO <- function(annot, gene_universe, effect){
  ## Create subuniverse:
  sub_universe <- gene_universe %>%
    subset(gene_universe$Name %in% unlist(annot$Parent))
  
  ## Run conditional hypergeometric test:
  runTestHypGeom <- function(sub_universe, onto){
    ## Constructing a GOHyperGParams objects or KEGGHyperGParams objects from a GeneSetCollection
    ## Then run hypergeometric test:
    GO_NO_fdr <- hyperGTest(GSEAGOHyperGParams(name="GO_set",
                                               geneSetCollection = gsc_universe,
                                               geneIds = as.vector(unique(sub_universe[["Name"]])), # gene ids for the selected gene set
                                               universeGeneIds = unique(gene_universe$Name),
                                               ontology = onto, # A string for GO to use ("BP", "CC", or "MF")
                                               pvalueCutoff = 0.05,
                                               conditional = TRUE, # see note above
                                               testDirection = "over")) # over represented GO terms
    
    # Use GOEnrich as a wrapper around GOStat for extra FDR comparison
    ## Does not solve all issues, but better than nothing. See: https://support.bioconductor.org/p/5571/
    GO_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                               gene.ids = as.vector(unique(sub_universe[["Name"]])),# genes in selected gene set
                                               univ.gene.ids = unique(gene_universe$Name),
                                               ontologies = onto, # A string for GO to use ("BP", "CC", or "MF")
                                               pvalue.cutoff = 0.05,
                                               cond = TRUE, # see note above
                                               test.dir = "over"),# over represented GO terms
                                  p.adjust.method = "fdr")
    
    
    return(list(GO_NO_fdr=GO_NO_fdr, GO_fdr=GO_fdr))
  }
  
  GO_MF <- runTestHypGeom(sub_universe = sub_universe, onto = "MF")
  GO_CC <- runTestHypGeom(sub_universe = sub_universe, onto = "CC")
  GO_BP <- runTestHypGeom(sub_universe = sub_universe, onto = "BP")
  
  # Get percentage of genes over reppresented in universe
  dfMFperc = GO_MF$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfCCperc = GO_CC$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfBPperc = GO_BP$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  
  # Add this information to FDR corrected table
  GO_MF_all = merge(GO_MF$GO_fdr, dfMFperc)
  GO_CC_all = merge(GO_CC$GO_fdr, dfCCperc)
  GO_BP_all = merge(GO_BP$GO_fdr, dfBPperc)
  
  # Merge the df MP and BP
  dfGO = rbind(GO_MF_all, GO_CC_all, GO_BP_all)
  
  dfGO = dfGO %>% mutate(Term = factor(x = GO.term, levels = GO.term),
                         Effect = factor(x = effect, levels = effect))
  
  # Relabel GO group names
  dfGO$GO.category[dfGO$GO.category %in% "CC"]="Cellular components"
  dfGO$GO.category[dfGO$GO.category %in% "BP"]="Biological processes"
  dfGO$GO.category[dfGO$GO.category %in% "MF"]="Molecular functions"
  
  return(dfGO)
}

##############
## Get DMS and subunite functions by BP

## Subselect those DMS present in at least 4 out of 8 BP, return a LIST with BP as names
get2keep = function(Compa, NBP = 4){
  x <- lapply(myPosList, function(x){unlist(x[[paste0("DMS_15pc_BP_", Compa)]])})
  f <- table(unlist((x))) # each DMS present between 1 and 8 times
  tokeep <- names(f)[f >= NBP]
  # print(length(tokeep))
  ## Keep the DMS present in 4 families minimum
  DMSBPlist_INTER4 <- lapply(x, function(x){x[x %in% tokeep]})
  ## Reorder by family:
  DMSBPlist_INTER4 <- DMSBPlist_INTER4[names(DMSBPlist_INTER4)[order(names(DMSBPlist_INTER4))]]
  return(DMSBPlist_INTER4)
}

# return a LIST with 1. dataframe of differential methylation and 2.subset of unite methylation raw object
get_dms.diffmeth.per1compa_4BPmin <- function(Compa){
  ## Extract a data frame with DMS and meth.diff for one given comparison:
  df.dms.methdiff = lapply(
    lapply(DMSBPlist, lapply, function(x){data.frame(DMS = paste(x$chr, x$end), meth.diff = x$meth.diff)}), 
    function(x){x[[paste0("DMS_15pc_BP_", Compa)]]})
  
  ## Add BP in the name of the column containing meth.diff:
  for (i in 1:length(names(df.dms.methdiff))){
    names(df.dms.methdiff[[i]])[names(df.dms.methdiff[[i]]) %in% "meth.diff"] = paste0("meth.diff_", names(df.dms.methdiff)[i])
  }
  
  # merge all elements of the list in one big data frame
  df.dms.methdiff = Reduce(function(...) merge(..., all=T), df.dms.methdiff)
  
  # keep DMS found in at least 4 BP:
  df.dms.methdiff = df.dms.methdiff[rowSums(!is.na(df.dms.methdiff[2:ncol(df.dms.methdiff)]))>= 4, ]
  
  # Subselect the original unite object for these DMS only 
  subUnite = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, 
                               which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, uniteCov14_G2_woSexAndUnknowChrOVERLAP$end) %in% df.dms.methdiff$DMS))
  return(list(df.dms.methdiff=df.dms.methdiff, subUnite=subUnite))
}

# sanity test that both functions give same results
# for (i in 1:4){
#   a=unique(unlist(get2keep(Compa = vecCompa[i])))
#   b=get_dms.diffmeth.per1compa_4BPmin(Compa = vecCompa[i])$df.dms.methdiff$DMS
#   print(table(a%in%b)); print(table(b%in%a))
# }
# 
# makePlotsobservedReactionNorms <- function(){
#   ### Extract differential methylation per comparison AND raw methylation values
#   ## G1 effect
#   A = get_dms.diffmeth.per1compa_4BPmin(Compa = vecCompa[1]) 
#   B = get_dms.diffmeth.per1compa_4BPmin(Compa = vecCompa[2])
#   ## G2 effect
#   C = get_dms.diffmeth.per1compa_4BPmin(Compa = vecCompa[3])
#   D = get_dms.diffmeth.per1compa_4BPmin(Compa = vecCompa[4])
#   
#   # Subselect the original unite object for all our DMS of interest
#   DMSofInterest = unique(c(A$df.dms.methdiff$DMS, B$df.dms.methdiff$DMS, C$df.dms.methdiff$DMS, D$df.dms.methdiff$DMS))
#   subUniteofInterest = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, 
#                                          which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, 
#                                                      uniteCov14_G2_woSexAndUnknowChrOVERLAP$end) %in% DMSofInterest))
#   
#   getObservedReacNorm <- function(DF, mytitle){
#     # Extract the top 5 more differentially methylated sites in this group
#     top5 = DF[apply(DF[2:ncol(DF)],1, mean, na.rm=T) %>% abs() %>% order(decreasing = T) %>% head(5),]
#     
#     # Get raw methylation values at these positions:
#     meth=methylKit::select(subUniteofInterest, which(paste(subUniteofInterest$chr, subUniteofInterest$end) %in% top5$DMS))
#     
#     if (nrow(meth) !=5){
#       print("ERROR!! Some top DMS found in several comparisons")
#     }
#     
#     dfmeth = meth%>% percMethylation()%>% data.frame()
#     dfmeth$DMS = paste(meth$chr, meth$end)
#     dfmeth=melt(dfmeth)
#     dfmeth$SampleID = as.character(dfmeth$variable)
#     
#     # Add brother pair and treatment info
#     dfmeth=merge(dfmeth, fullMetadata_OFFS[c("SampleID", "brotherPairID", "outcome", "patTrt")])
#     
#     # Make reaction norms plot per brother pair (expected: flat)
#     mean_data <- dfmeth %>% group_by(patTrt, outcome, brotherPairID, DMS) %>%
#       dplyr::summarize(value = mean(value, na.rm = TRUE)) %>% 
#       # Add bands of grey per chromosome for plot:
#       mutate(type=ifelse(as.numeric(as.factor(DMS))%%2, "A" , "B")) %>%  
#       data.frame()
#     
#     ggplot(mean_data, aes(x=outcome, y=value))+
#       facet_grid(DMS~brotherPairID) +
#       geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
#       scale_fill_manual(values = c("white", "grey"),guide="none")+
#       geom_point(aes(colour=patTrt))+
#       geom_line(aes(group=patTrt, colour=patTrt))+
#       scale_color_manual("Paternal (G1) treatment", values = c("black", "red"))+
#       xlab("Offspring (G2) treatment")+
#       ylab("Methylation value")+
#       ggtitle(mytitle)
#   }
#   
#   ##################
#   ## 1. Paternal effect only
#   G1dfonly = rbind(A$df.dms.methdiff[A$df.dms.methdiff$DMS %in% caseVennG1only,],
#                    B$df.dms.methdiff[B$df.dms.methdiff$DMS %in% caseVennG1only,])
#   plot1 = getObservedReacNorm(G1dfonly, mytitle = "Paternal effect only")
#   
#   ## 2. Offspring effect only
#   G2dfonly = rbind(C$df.dms.methdiff[C$df.dms.methdiff$DMS %in% caseVennG2only,],
#                    D$df.dms.methdiff[D$df.dms.methdiff$DMS %in% caseVennG2only,])
#   plot2 = getObservedReacNorm(G2dfonly, mytitle = "Offspring effect only")
#   
#   ## 3. G1 G2 NO interactions
#   G1G2NOinter_df = rbind(C$df.dms.methdiff[C$df.dms.methdiff$DMS %in% caseVennG1G2NOinter,],
#                          D$df.dms.methdiff[D$df.dms.methdiff$DMS %in% caseVennG1G2NOinter,])
#   plot3 = getObservedReacNorm(G1G2NOinter_df, mytitle = "G1 + G2 effect only")
#   
#   ## 4. G1 G2 WITH interactions
#   G1G2inter_df = rbind(C$df.dms.methdiff[C$df.dms.methdiff$DMS %in% caseVennG1G2inter,],
#                        D$df.dms.methdiff[D$df.dms.methdiff$DMS %in% caseVennG1G2inter,])
#   plot4 = getObservedReacNorm(G1G2inter_df, mytitle = "G1 : G2 effect only")
#   
#   ## 5. caseVennG2interNOG1s
#   G2interNOG1_df = rbind(C$df.dms.methdiff[C$df.dms.methdiff$DMS %in% caseVennG2interNOG1,],
#                          D$df.dms.methdiff[D$df.dms.methdiff$DMS %in% caseVennG2interNOG1,])
#   plot5 = getObservedReacNorm(G2interNOG1_df, mytitle = "G2 inter NO G1 effect only")
#   return(list(plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4, plot5=plot5))
# }

###########################
# Get association between DMS and gene
getAnnotDMS <- function(DMSvec, pasted=FALSE){
  if (pasted == FALSE){
    DMSvec =  paste(DMSvec$chr, DMSvec$end)
  }
  # Change the vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=sapply(strsplit(DMSvec, " "), `[`, 1), 
                                                  start=sapply(strsplit(DMSvec, " "), `[`, 2),
                                                  end=sapply(strsplit(DMSvec, " "), `[`, 2)))
  # nbrBP=DMSdf[[2]]), keep.extra.columns = T) # add brother pairs number
  A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = annotBed12)
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((A@dist.to.TSS$dist.to.feature>10000 | A@dist.to.TSS$dist.to.feature< -10000) &
                    rowSums(A@members) %in% 0)
  if (is_empty(rows2rm)){
    GRangeOBJ = GRangeOBJ
  } else {
    GRangeOBJ = GRangeOBJ[-rows2rm,]
  }
  ## Re-annotate the subsetted object
  B = annotateWithGeneParts(as(GRangeOBJ,"GRanges"),annotBed12)
  ## Get genes associated with these
  C = getAssociationWithTSS(B)
  # In which gene is the DMS?
  dmsdf = data.frame(DMS = paste(GRangeOBJ@seqnames, GRangeOBJ@ranges),
                     feature.name = C$feature.name)
  ## Get annotations for these genes
  subAnnot = data.frame(subset(annotGff3, Name %in% C$feature.name))
  
  # Extract gene symbol from the "Note" attribute
  subAnnot$Note = unlist(subAnnot$Note)
  subAnnot$GeneSymbol = str_extract(subAnnot$Note, "(?<=Similar to )(\\w+)")
  
  ## Associate DMS and gene
  featuredf = data.frame(feature.name = subAnnot$ID,
                         Note = subAnnot$Note,
                         GeneSymbol = subAnnot$GeneSymbol)
  
  return(merge(dmsdf, featuredf))
}

# ## Plot the most important CpG associated with significant axis 
# plotCpGexample_PCAsignif = function(DMSvec, mygene){
#   
#   # Get association between DMS and gene
#   DMSvsGene = getAnnotDMS(DMSvec)
#   
#   # Annotate the DMS
#   annotPCA <- getAnnotationFun(DMSdf = paste(DMSvec$chr, DMSvec$end), annotBed12 = annotBed12,
#                                annotGff3 = annotGff3, isDMDaDataframeWithBP = FALSE)
#   
#   # Check that all sequences with no name are "Protein of unknown function"
#   table(annotPCA$Note[is.na(annotPCA$GeneSymbol)] %in% "Protein of unknown function") ## all true if ok
#   
#   # Find DMS associated with a gene linked with PCA axis
#   dfDMS1gene = merge(annotPCA[annotPCA$GeneSymbol %in% mygene,],DMSvsGene) 
#   
#   raw=methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP,
#                         which(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr %in% dfDMS1gene$chrom & 
#                                 uniteCov14_G2_woSexAndUnknowChrOVERLAP$start %in% sapply(strsplit(dfDMS1gene$DMS," "), `[`, 2)))
#   
#   dfPlot = data.frame(chr=raw$chr, pos=raw$end)
#   dfPlot = cbind(dfPlot, data.frame(percMethylation(raw)))
#   dfPlot = melt(dfPlot, id.vars = c("chr", "pos")) %>% dplyr::rename("SampleID" = "variable")
#   
#   # Add sample group
#   dfPlot = merge(dfPlot, fullMetadata_OFFS[c("SampleID", "trtG1G2", "outcome", "patTrt", "brotherPairID")])
#   
#   # Plot
#   dfPlotSum = dfPlot %>% group_by(pos, outcome, patTrt, brotherPairID) %>%
#     dplyr::summarise(meanMeth=mean(value, na.rm=T)) %>% data.frame()
#   
#   ## All BP
#   ggplot(dfPlotSum)+
#     geom_point(aes(x=outcome, y=meanMeth, col=patTrt, group=patTrt), size = 3) +
#     geom_line(aes(x=outcome, y=meanMeth, col=patTrt, group=patTrt))+
#     scale_color_manual("Paternal treatment", values = c("black", "red"))+
#     # all data:
#     geom_point(data = dfPlot, aes(x=outcome, y=value, col=patTrt, group=patTrt), alpha=.5) +
#     facet_grid(pos~brotherPairID) +
#     ggtitle(paste("DMS in gene", annotPCA[annotPCA$GeneSymbol %in% mygene,"GeneSymbol"], ":", annotPCA[annotPCA$GeneSymbol %in% mygene,"description"]))+
#     xlab("Offspring treatment")+
#     ylab("Mean methylation")
# }
