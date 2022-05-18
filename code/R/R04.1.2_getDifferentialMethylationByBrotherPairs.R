## Differential methylation analyses BY BROTHER PAIR
## A. Balard
## February 2022

machine="mythinkpad" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
source("R02.3_DATALOAD.R")

######### 
## Calculate parental DMS (parDMS) for each brother pair

## CpG covered in half trt group, parental DMS:
nrow(uniteCov6_G1_woSexAndUnknowChrOVERLAP) # methylBase object with 1001880 rows

## Calculate DMS/DMR


#############
# Select only G1 fish from BP which have offspring
metadata_G1 <- fullMetadata_PAR[
  fullMetadata_PAR$brotherPairID %in% fullMetadata_OFFS$brotherPairID]

#### We will apply the following function to all BP:
vecBP <- unique(metadata_G1$brotherPairID)
# some parental BP didn't have offspring with sequencing, select them out
vecBP <- vecBP[vecBP %in% fullMetadata_OFFS$brotherPairID]

getG2atParDMSperBP <- function(BP){
  metadata1BP_G1 <- metadata_G1[metadata_G1$brotherPairID %in% BP,]
  
  # Select methylation object for 1 BP from methylobject with CpG covered in ALL fish
  uniteCov1BP_G1 <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                               treatment = metadata1BP_G1$trtG1G2_NUM,
                               sample.ids = metadata1BP_G1$ID)
  
  DMS15pc_1BP <- getDiffMethSimple(uniteCov1BP_G1, metadata1BP_G1)
  ## NOTE: performing 'fast.fisher' instead of 'F' for two groups testing.
  
  #############
  ## Get methylation in offspring of these 2 brothers at the parDMS positions:
  
  # metadata
  metadata1BP_G2 <- fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP,]
  
  # methyldata
  uniteCov1BP_G2 <- reorganize(methylObj = uniteCovALL_woSexAndUnknowChr,
                               treatment = metadata1BP_G2$trtG1G2_NUM,
                               sample.ids = metadata1BP_G2$ID)
  
  # keep only CpG corresponding to parDMS in this BP, split in hypo/hyper
  parDMS1BP_G2_hypo <- methylKit::select(uniteCov1BP_G2, 
                                         which(paste(uniteCov1BP_G2$chr, uniteCov1BP_G2$start) %in% 
                                                 paste(DMS15pc_1BP[DMS15pc_1BP$meth.diff<0, ]$chr, 
                                                       DMS15pc_1BP[DMS15pc_1BP$meth.diff<0, ]$start)))
  parDMS1BP_G2_hyper <- methylKit::select(uniteCov1BP_G2, 
                                          which(paste(uniteCov1BP_G2$chr, uniteCov1BP_G2$start) %in% 
                                                  paste(DMS15pc_1BP[DMS15pc_1BP$meth.diff>0, ]$chr, 
                                                        DMS15pc_1BP[DMS15pc_1BP$meth.diff>0, ]$start)))
  
  # Get average methylation per individual at these CpGs
  aveMeth_hypo <- apply(percMethylation(parDMS1BP_G2_hypo), 2, mean)
  aveMeth_hyper <- apply(percMethylation(parDMS1BP_G2_hyper), 2, mean)
  
  # make datasets for statistical test
  dfhypoG2parDMS <- data.frame(SampleID=names(aveMeth_hypo), aveMeth=aveMeth_hypo)
  dfhyperG2parDMS <- data.frame(SampleID=names(aveMeth_hyper), aveMeth=aveMeth_hyper)
  
  # add metadata
  dfhypoG2parDMS <- merge(dfhypoG2parDMS, fullMetadata_OFFS)
  dfhyperG2parDMS <- merge(dfhyperG2parDMS, fullMetadata_OFFS)
  
  # reorder sample ID cause it was mixed up
  dfhypoG2parDMS = dfhypoG2parDMS[order(as.numeric(gsub("S", "", dfhypoG2parDMS$ID))),]
  dfhyperG2parDMS = dfhyperG2parDMS[order(as.numeric(gsub("S", "", dfhyperG2parDMS$ID))),]
  
  return(list(parDMS1BP_G2_hypo=parDMS1BP_G2_hypo, parDMS1BP_G2_hyper=parDMS1BP_G2_hyper,
              dfhypoG2parDMS=dfhypoG2parDMS, dfhyperG2parDMS=dfhyperG2parDMS))
} 

## Here all the dataset for each BP are stored (metadata for univariate stats/methylKit object for multivariate)
results <- lapply(vecBP, getG2atParDMSperBP)

######### 
## 1. Multivariate analysis: clustering of samples by trt
results4PCA <- results # we modify the object
for (i in 1:8){
  results4PCA[[i]]$parDMS1BP_G2_hypo@treatment = 
    merge(data.frame(old = results4PCA[[i]]$parDMS1BP_G2_hypo@treatment),
          data.frame(old=c(2,3,5,6), new=1:4))$new
  
  # change number for trt to visualise
  results4PCA[[i]]$parDMS1BP_G2_hypo@sample.ids = results4PCA[[i]]$dfhypoG2parDMS$Tr
  
  PCASamples(results4PCA[[i]]$parDMS1BP_G2_hypo)
  Sys.sleep(10)
}

######### 
## 2. Univariate analysis: lmer based on average methylation level at parDMS in each BP
aveMethStats <- function(BPnum, hypoOrhyper){
  if (hypoOrhyper=="hypo"){
    df = results[[BPnum]]$dfhypoG2parDMS
  } else if (hypoOrhyper=="hyper"){
    df = results[[BPnum]]$dfhyperG2parDMS
  }
  
  modfull <- lmer(aveMeth ~ patTrt * outcome + (1|Sex), data = df)
  modnointer <- lmer(aveMeth ~ patTrt + outcome + (1|Sex), data = df)
  
  # test significance of factors
  lrt_inter = lrtest(modfull, modnointer)
  isInterSignif = lrt_inter$`Pr(>Chisq)`[2]<0.05
  
  ## Prediction plot
  p1 = plot_model(modfull, type = "pred", terms = c("patTrt", "outcome"))+
    theme_cleveland()+
    ggtitle(label = "Predicted average methylation", 
            subtitle = paste0("parental DMS ", hypoOrhyper, "methylated\nupon infection (", vecBP[BPnum], ")"))+
    xlab("Paternal treatment")+
    labs(col="Offspring trt: ")+
    scale_color_manual(values = c("black", "red"))
  if(isInterSignif){
    p1=p1 + theme(panel.background=element_rect(colour="yellow", size = 5))}
  
  ## Extract slope between CC-CT and TC-TT in each BP
  pred=data.frame(ggpredict(modfull,terms = c("patTrt", "outcome")))
  CC_CTslope = pred$predicted[pred$x %in% "controlP" & pred$group %in% "infected"] - 
    pred$predicted[pred$x %in% "controlP" & pred$group %in% "control"]
  TC_TTslope = pred$predicted[pred$x %in% "infectedP" & pred$group %in% "infected"] - 
    pred$predicted[pred$x %in% "infectedP" & pred$group %in% "control"]
  
  ## average in parent group
  avePred=data.frame(ggeffect(modfull,terms = c("patTrt")))
  CparMean=avePred$predicted[avePred$x %in% "controlP"]
  TparMean=avePred$predicted[avePred$x %in% "infectedP"]
  
  mysummary=data.frame(BP = unique(df$brotherPairID), 
                       isInterSignif=isInterSignif, hypoOrhyper=hypoOrhyper,
                       CC_CTslope=CC_CTslope, TC_TTslope=TC_TTslope,
                       CparMean = CparMean,TparMean=TparMean)
  
  return(list(p1=p1, mysummary=mysummary))
}

## Plot for all BP
plist <- list() # empty plot list
for (i in 1:length(vecBP)){
  plist[[i]] <- aveMethStats(BPnum = i, "hypo")$p1
} 
do.call(grid.arrange, c(plist, ncol = 4))

plist <- list() # empty plot list
for (i in 1:length(vecBP)){
  plist[[i]] <- aveMethStats(BPnum = i, "hyper")$p1
} 
do.call(grid.arrange, c(plist, ncol = 4))

## Get dataframe of summarised values for all BP
hypodf=do.call(rbind, lapply(1:length(vecBP), 
                             function(x) aveMethStats(BPnum = x, hypoOrhyper = "hypo")$mysummary))
hypodf=melt(hypodf, id.vars = c("BP", "isInterSignif", "hypoOrhyper"))

hyperdf=do.call(rbind, lapply(1:length(vecBP), 
                              function(x) aveMethStats(BPnum = x, hypoOrhyper = "hyper")$mysummary))
hyperdf=melt(hyperdf, id.vars = c("BP", "isInterSignif", "hypoOrhyper"))

## plot summarised values
## reaction norm = pure transgenerational effect
ggplot(hypodf[hypodf$variable %in% c("CC_CTslope", "TC_TTslope"),], 
       aes(x=variable, y=value))+
  geom_point(aes(col=isInterSignif), size = 3)+
  theme_cleveland()+
  scale_color_manual(values = c("black", "red"))+
  geom_line(aes(group=BP,col=isInterSignif))
ggplot(hyperdf[hyperdf$variable %in% c("CC_CTslope", "TC_TTslope"),], 
       aes(x=variable, y=value))+
  geom_point(aes(col=isInterSignif), size = 3)+
  theme_cleveland()+
  scale_color_manual(values = c("black", "red"))+
  geom_line(aes(group=BP,col=isInterSignif))

## marginal prediction of parental infection = transgenerational effect + genetic brothers
ggplot(hypodf[hypodf$variable %in% c("CparMean", "TparMean"),], 
       aes(x=variable, y=value))+
  geom_point(aes(col=isInterSignif), size = 3)+
  theme_cleveland()+
  scale_color_manual(values = c("black", "red"))+
  geom_line(aes(group=BP,col=isInterSignif))
ggplot(hyperdf[hyperdf$variable %in% c("CparMean", "TparMean"),], 
       aes(x=variable, y=value))+
  geom_point(aes(col=isInterSignif), size = 3)+
  theme_cleveland()+
  scale_color_manual(values = c("black", "red"))+
  geom_line(aes(group=BP,col=isInterSignif))



