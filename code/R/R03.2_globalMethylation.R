## Global methylation analyses
## A. Balard
## November 2021

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
#machine="apocrita"
machine="mythinkpad"
## Load methylation data
loadALL = TRUE # load all uniteCov objects
source("R02.2_loadMethyldata.R")

#############################################################
### PART 1: Methylation profiles, CpG present in all fish ###
#############################################################

## Dendogram of methylations
## All samples:
pdf("Rfigures/clusterALLCpG.pdf", width = 16, height = 7)
makePrettyMethCluster(uniteCovALL_woSexAndUnknowChr, fullMetadata,
                      my.cols.trt=c("#333333ff","#ff0000ff","#ffe680ff","#ff6600ff","#aaccffff","#aa00d4ff"),
                      my.cols.fam = c(1:4))
dev.off()

## offspring: (add TRUE to script loadMethylData)
pdf("Rfigures/clusterALLCpG_offspings.pdf", width = 17, height = 7)
makePrettyMethCluster(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS,
                      my.cols.trt=c("#ffe680ff","#ff6600ff", "#aaccffff", "#aa00d4ff"),
                      my.cols.fam = c(1:4))
dev.off()

uniteCovALL_G2_woSexAndUnknowChr@treatment
fullMetadata_OFFS$trtG1G2_NUM


############################
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

###################
## Let's run Adonis for offspring (considering CpG shared by ALL)
myadonisFUN(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
# PAT        1 0.002782 0.01470 1.6344 0.000999 ***
#   outcome    1 0.001909 0.01009 1.1216 0.068931 .  
# Sex        1 0.002418 0.01277 1.4202 0.010989 *  

# trtG1G2    3 0.006418 0.03391 1.2568 0.000999 ***
#   Sex        1 0.002407 0.01272 1.4141 0.008991 ** 

########## NMDS
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

#### RUN Goodness of fit
# myGOF.NMDS.FUN(dataset = uniteCovALL_G2_woSexAndUnknowChr) # Goodness of fit for NMDS 
# suggested the presence of six dimensions with a stress value <0.1 and 2 with > 0.2

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

NMDSanalysis <- myNMDS(dataset = uniteCovALL_G2_woSexAndUnknowChr, metadata = fullMetadata_OFFS)
NMDSanalysis$NMDSplot
save(NMDSanalysis, file = "../../data/fig/NMDSplots.RData")
