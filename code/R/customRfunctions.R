################################################
## List of custom functions used in the project:
## Some from Melanie heckwolf & Charley (Eugenie) Yen

# deroman() turns roman numbers in arabic ones in genome chr names
# makePrettyMethCluster() to make a cluster figure with annotation bars

## Adonis functions:
# makePercentMetMat() creates a matrix containing percent methylation values
# makeDatadistFUN() to make a pairwise distance matrix for ADONIS
# AdonisWithinG1trtFUN() to run adonis within G1 treatment (father control or infected)

## NMDS functions:
# myGOF.NMDS.FUN check goodness of fit
# myNMDSFUN runs my NMDS

## GO functions:
# makedfGO() takes an annotation df, a gene universe and the name of an effect
# and returns a GO df
# makeGOplot() makes a GO plot from dfGO
# makeGOplotslim() makes a GO plot with slim GO terms from dfGO
# getGo2Goslim() map a GOslim term to its GO terms

## PCA functions:
# getG2methAtCpGs() Get methylation at list of CpG for all offspring
# getResPCA() run PCA with imputation of missing data
# getModPCA() test model from PCA
# plotScreePlot() get scree plot
# myPlotPCA() plot PCA
# plotModBCI_PCA() plot the model

##############

#the stickleback genome is annotated as chrI, chrII, chrIII ... to chrUn (unknown, includes all scaffolds)
#this function takes the roman rumbers and turns them into actual numbers:
deroman <- function(x){ x %>% str_remove(.,"Gy_chr") %>%
    ifelse(. %in% c("Un","M"),., as.roman(.) %>%
             as.integer() %>% as.character() %>% str_pad(.,2,"left",0))
}

makePrettyMethCluster <- function(OBJ, metadata, my.cols.trt, my.cols.fam, nbrk, rect = T){
  ## Reorder metadata by sample ID, as OBJ methylkit!
  metadata = metadata[order(as.numeric(gsub("S", "", metadata$SampleID))),]
  
  ## Check
  if (!is.na(table(OBJ@sample.ids == metadata$SampleID)["FALSE"])){
    stop("check the samples order or similarity before both methylkit and metadata objects!")
  }
  
  ## To add color bars:
  # Generate color palette
  pal = qualpalr::qualpal(55, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))
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
  
  ## Make dendrogram
  mydendro <- clusterSamples(OBJ, dist="correlation", method="ward", plot=FALSE)
  dend = as.dendrogram(mydendro)
  
  ## and plot
  dend %>% plot(main=paste(OBJ@context, "methylation clustering\n",
                           "Distance method: correlation; Clustering method: ward.D"),
                font.main = 1, cex.main = .8, ylab = "Height", nodePar = list(lab.cex = 0.6, pch = c(NA, NA)))
  if(rect){
    dend %>% rect.dendrogram(k=nbrk, border = 8, lty = 5, lwd = 2)  
  }
  colored_bars(cbind(col_trt, col_trtPAT, brotherPairID), dend, y_shift = -0.1, #col_clutch
               rowLabels = c("G2 treatment", "G1 treatment", "G1 family")) #"Clutch"
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
  if (length((which(SD<0.3)))!=0){
    perc.meth = perc.meth[-which(SD<0.3),]  
  }
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

########
## GO ##
########
makedfGO <- function(annot, gene_universe, effect, label){
  ## Create subuniverse:
  sub_universe <- gene_universe %>%
    subset(gene_universe$Name %in% unlist(annot$Parent))
  
  ## Run conditional hypergeometric test:
  runTestHypGeom <- function(sub_universe, onto){
    ## Constructing a GOHyperGParams objects or KEGGHyperGParams objects from a GeneSetCollection
    ## Then run hypergeometric test:
    GO_NO_fdr <- hyperGTest(GSEAGOHyperGParams(
      name="GO_set",
      geneSetCollection = gsc_universe,
      geneIds = as.vector(unique(sub_universe[["Name"]])), # gene ids for the selected gene set
      universeGeneIds = unique(gene_universe$Name),
      ontology = onto, # A string for GO to use ("BP", "CC", or "MF")
      pvalueCutoff = 0.05,
      conditional = TRUE, # see note above
      testDirection = "over")) # over represented GO terms
    
    # Use GOEnrich as a wrapper around GOStat for extra FDR comparison
    ## Does not solve all issues, but better than nothing. See: https://support.bioconductor.org/p/5571/
    GO_fdr <- joinGOEnrichResults(goEnrichTest(
      gsc=gsc_universe,
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
  
  # Get percentage of genes over represented in universe
  dfMFperc = GO_MF$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>%
    dplyr::select(c("Term", "genePercent","Count", "Size")) %>% dplyr::rename(GO.name=Term)
  dfCCperc = GO_CC$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>%
    dplyr::select(c("Term", "genePercent","Count", "Size")) %>% dplyr::rename(GO.name=Term)
  dfBPperc = GO_BP$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>%
    dplyr::select(c("Term", "genePercent","Count", "Size")) %>% dplyr::rename(GO.name=Term)
  
  # Add this information to FDR corrected table
  GO_MF_all = merge(GO_MF$GO_fdr, dfMFperc)
  GO_CC_all = merge(GO_CC$GO_fdr, dfCCperc)
  GO_BP_all = merge(GO_BP$GO_fdr, dfBPperc)
  
  # Merge the df MP and BP
  dfGO = rbind(GO_MF_all, GO_CC_all, GO_BP_all)
  
  dfGO = dfGO %>% mutate(Term = factor(x = GO.term, levels = GO.term),
                         Label = factor(x = label, levels = label),
                         Effect = factor(x = effect, levels = effect))
  
  # Relabel GO group names
  dfGO$GO.category[dfGO$GO.category %in% "CC"]="Cellular components"
  dfGO$GO.category[dfGO$GO.category %in% "BP"]="Biological processes"
  dfGO$GO.category[dfGO$GO.category %in% "MF"]="Molecular functions"
  
  return(dfGO)
}

## make a GO plot from dfGO
makeGOplot <- function(dfGO, posleg="top"){
  dfGO %>%
    dplyr::filter(p.value.adjusted < 0.05) %>% 
    ggplot(aes(x=Effect, y = factor(GO.name))) +
    geom_point(aes(color = p.value.adjusted, size = genePercent)) +
    scale_color_gradient(
      name="adjusted\np-value", low = "red", high = "blue", 
      limits = c(0, 0.05), breaks = c(0, 0.02, 0.04), labels =c("0", "0.02", "0.04")) +
    scale_size_continuous(name = "% of genes")+
    theme_bw() + ylab("") + xlab("") +
    theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
          legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
          legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), # grey box for legend
          legend.position="top",
          axis.text.y = element_text(size = 8),  # Decrease y-axis text size
          axis.text.x = element_text(size = 8, hjust = 1)  # Increase x-axis text size 
    )+
    facet_wrap(.~fct_inorder(GO.category), scales = "free")
}

###############
## Check GO slim terms for easy interpretation
## GO subsets (also known as GO slims) are condensed versions of the GO containing a subset of the terms.
# dl the GO slim Developed by GO Consortium for the Alliance of Genomes Resources
# download.file(url = "https://current.geneontology.org/ontology/subsets/goslim_agr.obo",
# destfile = "../../data/goslim_agr.obo")
slim <- GSEABase::getOBOCollection("../../data/goslim_agr.obo")

## make a GO plot with slim GO terms from dfGO
makeGOslim <- function(dfGO, posleg="top"){
  
  myfun <- function(myeffect, myGOcat){
    # Create the GOCollection
    go_collection <- GOCollection(dfGO[
      dfGO$p.value.adjusted < 0.05 & dfGO$Effect %in% myeffect,"GO.term"])
    
    # Perform the GO slim mapping
    slimdf <- GSEABase::goSlim(idSrc = go_collection, 
                               slimCollection = slim,
                               ontology = myGOcat)
    
    # Map the original GO terms to the slim terms
    mappedIds <- function(df, collection, OFFSPRING) {
      map <- as.list(OFFSPRING[rownames(df)])
      mapped <- lapply(map, intersect, ids(collection))
      df[["mapped_go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L))
      
      # Get GO term names
      go_names <- AnnotationDbi::select(GO.db, keys = unlist(mapped), columns = "TERM", keytype = "GOID")
      go_names_list <- split(go_names$TERM, go_names$GOID)
      df[["mapped_go_names"]] <- vapply(mapped, function(x) {
        paste(go_names_list[x], collapse = ";")
      }, character(1L))
      
      # Get GO slim term full names
      slim_names <- AnnotationDbi::select(GO.db, keys = rownames(df), columns = "TERM", keytype = "GOID")
      df[["go_slim_full_name"]] <- slim_names$TERM
      
      df
    }
    
    # Use the appropriate OFFSPRING database based on the ontology
    offspring_db <- switch(myGOcat,
                           "BP" = GO.db::GOBPOFFSPRING,
                           "CC" = GO.db::GOCCOFFSPRING,
                           "MF" = GO.db::GOMFOFFSPRING)
    
    slimdf_with_terms <- mappedIds(slimdf, go_collection, offspring_db)
    
    slimdf_with_terms %>%
      filter(Count != 0) %>%
      dplyr::mutate(GO.category = myGOcat, Effect = myeffect)
  }
  
  list = list(NULL) ; i=1
  for (myeffect in unique(dfGO$Effect)){
    for (myGOcat in c("BP", "CC", "MF")){
      list[[i]] = myfun(myeffect, myGOcat)
      i = i+1
    }
  }
  
  dfGOslim = do.call(rbind, list)
  
  ## Order by percent
  dfGOslim = dfGOslim[order(dfGOslim$Percent), ]
  
  return(dfGOslim)
}

makeGOslimPlot <- function(dfGOslim, posleg="top"){
  
  dfGOslim %>%
    ggplot(aes(x=Effect, y = Term)) +
    geom_point(aes(size = Percent)) +
    scale_size_continuous(name = "% of GO terms in this GO slim category", range = c(1,8))+
    theme_bw() + ylab("") + xlab("") +
    theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
          legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
          legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), # grey box for legend
          legend.position=posleg,
          axis.text.y = element_text(size = 8),  # Decrease y-axis text size
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1)  # Increase x-axis text size and rotate
    )+
    facet_grid(Effect~fct_inorder(GO.category), scales="free",space = "free")+
    coord_flip() 
}

#########
## PCA ##
#########
## Get methylation at CpG for offspring
getG2methAtCpGs <- function(DMSlist){
  obj=methylKit::select(
    uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, 
    which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, 
                uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start, sep = "_") %in% 
            DMSlist))
  
  DMS=paste(obj$chr, obj$start, sep = "_")
  
  obj=percMethylation(obj) %>% data.frame 
  obj$DMS = DMS
  
  obj=melt(obj)
  names(obj)[names(obj)%in% "variable"]="SampleID"
  names(obj)[names(obj)%in% "value"]="G2methylation"
  obj=merge(obj, fullMetadata_OFFS[c("SampleID", "brotherPairID", "Sex",
                                     "trtG1G2", "patTrt", "No.Worms", "BCI")])
  obj$patTrt[obj$patTrt %in% "controlP"]="Control"
  obj$patTrt[obj$patTrt %in% "infectedP"]="Exposed"
  
  ## Pivot to have one column per DMS
  obj = obj %>% pivot_wider(names_from = DMS, values_from = G2methylation) %>%
    data.frame()
  
  ## reorder by sample original level
  obj = obj[order(obj$SampleID, levels(obj$SampleID)),]
  
  return(obj)
}

getResPCA <- function(x, incomplete = T){
  
  ## make samples as row names
  x = x %>% `rownames<-`(.[,1]) %>% 
    dplyr::select(matches("Gy_chr"))
  ## check samples order
  print("Samples order:") ; print(head(rownames(x)))
  
  if (incomplete==TRUE){
    # estimate the number of components from incomplete data
    nb <- estim_ncpPCA(x, scale = T)
    # impute the table
    res.comp <- imputePCA(x, ncp = nb$ncp, scale = T)
    x = res.comp$completeObs
  }
  # 2. run PCA
  res.PCA = FactoMineR::PCA(x, scale.unit = T, graph = FALSE) # perform PCA
  return(res.PCA)
}

getModPCA <- function(res.PCA){
  # extract axes 1, 2
  dfPCA = data.frame(PCA1=res.PCA$ind$coord[,1], # axis 1
                     PCA2=res.PCA$ind$coord[,2], # axis 2
                     SampleID = row.names(res.PCA$ind$coord))
  # add metadata
  dfPCA = merge(dfPCA, fullMetadata_OFFS)                     
  ## reorder samples following PCA
  dfPCA = dfPCA[match(rownames(res.PCA$ind$coord), dfPCA$SampleID),]
  
  ## check samples order
  print("Samples order:") ; print(head(dfPCA$SampleID))
  
  # Correlation with parasite load/BCI
  mod = lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), 
             data=dfPCA)
  
  print(lmerTest::step(mod, reduce.random = F))
  
  ## Model selection:
  modSel = lmer(formula = attr(attr(lmerTest::step(mod, reduce.random = F), "drop1"), "heading")[3],
                data=dfPCA, REML = F)
  
  return(list(modSel = modSel, metadata = dfPCA))
}

getModPCA_bygroup <- function(res.PCA){
  # extract axes 1, 2
  dfPCA = data.frame(PCA1=res.PCA$ind$coord[,1], # axis 1
                     PCA2=res.PCA$ind$coord[,2], # axis 2
                     SampleID = row.names(res.PCA$ind$coord))
  # add metadata
  dfPCA = merge(dfPCA, fullMetadata_OFFS)                     
  ## reorder samples following PCA
  dfPCA = dfPCA[match(rownames(res.PCA$ind$coord), dfPCA$SampleID),]
  
  ## check samples order
  print("Samples order:") ; print(head(dfPCA$SampleID))
  
  # Correlation with parasite load/BCI
  mod = lmer(BCI ~ PCA1*PCA2*No.Worms + (1|brotherPairID)+ (1|Sex), 
             data=dfPCA)
  
  print(lmerTest::step(mod, reduce.random = F))
  
  ## Model selection:
  modSel = lmer(formula = attr(attr(lmerTest::step(mod, reduce.random = F), "drop1"), "heading")[3],
                data=dfPCA, REML = F)
  
  return(list(modSel = modSel, metadata = dfPCA))
}

plotScreePlot <- function(resPCA){
  barplot(resPCA$eig[, 2], names.arg=1:nrow(resPCA.allDMS$eig), 
          main = "Variances",
          xlab = "Principal Components",
          ylab = "Percentage of variances",
          col ="steelblue")
  # Add connected line segments to the plot
  lines(x = 1:nrow(resPCA$eig), resPCA.allDMS$eig[, 2], 
        type="b", pch=19, col = "red")
}

myPlotPCA <- function(resPCA, resMod_PCA){
  fviz_pca_ind(resPCA, label="none", 
               habillage=resMod_PCA$metadata$trtG1G2,
               pointsize =3, addEllipses=TRUE)+
    scale_color_manual(values = colOffs)+
    scale_fill_manual(values = colOffs)
}

plotModBCI_PCA <- function(model, resMod_PCA, Signifterms){
  p1 = plot_model(modFULL)
  
  ## predicted and observed data
  pred_data <- ggpredict(modFULL, terms = Signifterms)  
  dfobs = resMod_PCA$metadata
  dfobs$x = dfobs[[Signifterms[1]]]
  dfobs$predicted = dfobs$BCI
  dfobs$group = dfobs$PAT
  dfobs$trtG1G2 = dfobs$trtG1G2
  
  p2 = ggplot(pred_data, aes(x = x, y = predicted)) +
    geom_line(aes(linetype = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
                alpha = 0.2, linetype = "blank", fill = "grey")+
    geom_jitter(data=dfobs, aes(fill = trtG1G2), size = 3, pch=21,col="black",
                height = 0, width = .1) +
    labs(x = Signifterms[1], y = "Body Condition Index", 
         color = "PAT", fill = "PAT")+
    scale_color_manual(values = setNames(colOffs, NULL)[1:2])+
    scale_fill_manual(values = colOffs)
  
  return(list(p1=p1, p2=p2))
}
