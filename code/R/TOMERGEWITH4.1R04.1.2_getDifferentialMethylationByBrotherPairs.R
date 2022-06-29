## Differential methylation analyses BY BROTHER PAIR
## A. Balard
## February 2022

machine="mythinkpad" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
loadannot = TRUE
sourceDMS = FALSE
source("R02.3_DATALOAD.R")

##########################
## I. Top-down approach ##
##########################

######### 
## Calculate parental DMS (parDMS) for each brother pair

## CpG covered in half trt group, parental DMS:
nrow(uniteCov6_G1_woSexAndUnknowChrOVERLAP) # methylBase object with 1001880 rows

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
  
  return(list(DMS15pc_1BP_G1 = DMS15pc_1BP,
              parDMS1BP_G2_hypo=parDMS1BP_G2_hypo, parDMS1BP_G2_hyper=parDMS1BP_G2_hyper,
              dfhypoG2parDMS=dfhypoG2parDMS, dfhyperG2parDMS=dfhyperG2parDMS))
} 

## Here all the dataset for each BP are stored (metadata for univariate stats/methylKit object for multivariate)
run=FALSE
if (run==TRUE){
  results <- lapply(vecBP, getG2atParDMSperBP)
}

######### 
## 1. Multivariate analysis: clustering of samples by trt
run=FALSE
if (run==TRUE){
  results4PCA <- results # we modify the object
  for (i in 1:8){
    results4PCA[[i]]$parDMS1BP_G2_hypo@treatment = 
      merge(data.frame(old = results4PCA[[i]]$parDMS1BP_G2_hypo@treatment),
            data.frame(old=c(2,3,5,6), new=1:4))$new
    
    # change number for trt to visualise
    results4PCA[[i]]$parDMS1BP_G2_hypo@sample.ids = results4PCA[[i]]$dfhypoG2parDMS$Tr
    ## Plot for all BP
    PCASamples(results4PCA[[i]]$parDMS1BP_G2_hypo)
    Sys.sleep(10)
  }
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
  CC_CTslope = pred[pred$x %in% "controlP" & pred$group %in% "infected", "predicted"] -
    pred[pred$x %in% "controlP" & pred$group %in% "control", "predicted"]
  TC_TTslope = pred[pred$x %in% "infectedP" & pred$group %in% "infected", "predicted"] -
    pred[pred$x %in% "infectedP" & pred$group %in% "control", "predicted"]
  
  ## average in parent group
  avePred=data.frame(ggeffects::ggpredict(modfull,terms = c("patTrt")))
  CparMean=avePred$predicted[avePred$x %in% "controlP"]
  TparMean=avePred$predicted[avePred$x %in% "infectedP"]
  
  mysummary=data.frame(BP = unique(df$brotherPairID),
                       isInterSignif=isInterSignif, hypoOrhyper=hypoOrhyper,
                       CC_CTslope=CC_CTslope, TC_TTslope=TC_TTslope,
                       CparMean = CparMean,TparMean=TparMean)
  
  return(list(p1=p1, mysummary=mysummary))
}

run=FALSE
if (run==TRUE){
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
  p1 = ggplot(hypodf[hypodf$variable %in% c("CC_CTslope", "TC_TTslope"),], 
              aes(x=variable, y=value))+
    geom_point(aes(col=isInterSignif), size = 3)+
    theme_cleveland()+
    scale_color_manual(values = c("black", "red"))+
    geom_line(aes(group=BP,col=isInterSignif))+
    geom_label_repel(data = hypodf[hypodf$variable %in% "TC_TTslope",], 
                     aes(x=variable, y=value, label=BP), xlim  = 2.1, segment.color="grey")
  p2 = ggplot(hyperdf[hyperdf$variable %in% c("CC_CTslope", "TC_TTslope"),], 
              aes(x=variable, y=value))+
    geom_point(aes(col=isInterSignif), size = 3)+
    theme_cleveland()+
    scale_color_manual(values = c("black", "red"))+
    geom_line(aes(group=BP,col=isInterSignif))+
    geom_label_repel(data = hyperdf[hyperdf$variable %in% "TC_TTslope",], 
                     aes(x=variable, y=value, label=BP), xlim  = 2.1, segment.color="grey")
  ggarrange(p1, p2, labels = c("hypo-parDMS", "hyper-parDMS"), ncol = 2)
  
  ## marginal prediction of parental infection = transgenerational effect + genetic brothers
  p1 = ggplot(hypodf[hypodf$variable %in% c("CparMean", "TparMean"),], 
              aes(x=variable, y=value))+
    geom_point(aes(col=isInterSignif), size = 3)+
    theme_cleveland()+
    scale_color_manual(values = c("black", "red"))+
    geom_line(aes(group=BP,col=isInterSignif))+
    geom_label_repel(data = hypodf[hypodf$variable %in% "TparMean",], 
                     aes(x=variable, y=value, label=BP), xlim  = 2.1, segment.color="grey")
  p2 = ggplot(hyperdf[hyperdf$variable %in% c("CparMean", "TparMean"),],
              aes(x=variable, y=value))+
    geom_point(aes(col=isInterSignif), size = 3)+
    theme_cleveland()+
    scale_color_manual(values = c("black", "red"))+
    geom_line(aes(group=BP,col=isInterSignif))+
    geom_label_repel(data = hyperdf[hyperdf$variable %in% "TparMean",], 
                     aes(x=variable, y=value, label=BP), xlim  = 2.1, segment.color="grey")
  ggarrange(p1, p2, labels = c("hypo-parDMS", "hyper-parDMS"), ncol = 2)
  
  ## Binomial test of our hypothesis: if the average methylation is HIGH in G2 from control G1,
  ## then the interaction G1trt:G2trt is not significant
  binom.test(x = 8, n = 8, p = 0.5, conf.level = 0.95) # p-value = 0.007812
  
  binom.test(x = 6, n = 8, p = 0.5, conf.level = 0.95) # p-value = 0.007812
} 

##############################################################
#### Compare per BP diffmeth in C-T, CC-CT and TC-TT groups ##
run=FALSE
if (run==TRUE){
  getDFpercentHyper <- function(i){
    res = results[[i]]
    
    ## Hypo
    PM=percMethylation(res$parDMS1BP_G2_hypo)
    ## Calculate average methylation for CC, CT, TC, TT
    trts = levels(fullMetadata_OFFS$trtG1G2)
    df=data.frame(
      A = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[1]]]),
      B = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[2]]]),
      C = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[3]]]),
      D = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[4]]]))
    names(df) = paste0("ave_", trts)
    ## Calculate difference Exposed minus Control in both parental group
    df = data.frame(DiffMeth_CT_NEG1 = df$ave_NE_exposed - df$ave_NE_control,
                    DiffMeth_CT_EG1 = df$ave_E_exposed - df$ave_E_control,
                    chr = res$parDMS1BP_G2_hypo$chr, start = res$parDMS1BP_G2_hypo$start, 
                    end = res$parDMS1BP_G2_hypo$end)
    ## Merge with parDMS values
    dfhypo = merge(methylKit::getData(res$DMS15pc_1BP_G1), df)
    
    ## Hyper
    PM=percMethylation(res$parDMS1BP_G2_hyper)
    ## Calculate average methylation for CC, CT, TC, TT
    trts = levels(fullMetadata_OFFS$trtG1G2)
    df=data.frame(
      A = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[1]]]),
      B = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[2]]]),
      C = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[3]]]),
      D = rowMeans(PM[,colnames(PM) %in% fullMetadata_OFFS$SampleID[
        fullMetadata_OFFS$trtG1G2 %in% trts[4]]]))
    names(df) = paste0("ave_", trts)
    ## Calculate difference Exposed minus Control in both parental group
    df = data.frame(DiffMeth_CT_NEG1 = df$ave_NE_exposed - df$ave_NE_control,
                    DiffMeth_CT_EG1 = df$ave_E_exposed - df$ave_E_control,
                    chr = res$parDMS1BP_G2_hyper$chr, start = res$parDMS1BP_G2_hyper$start, 
                    end = res$parDMS1BP_G2_hyper$end)
    ## Merge with parDMS values
    dfhyper = merge(methylKit::getData(res$DMS15pc_1BP_G1), df)
    
    dfhypohyper = rbind(dfhypo, dfhyper) # NB: positions NOT ordered
    
    result = data.frame(BP=unique(res$dfhypoG2parDMS$brotherPairID),
                        percentHyper_G1 = sum(dfhypohyper$meth.diff > 0) / nrow(dfhypohyper),
                        percentHyper_G2_NEG1 = sum(dfhypohyper$DiffMeth_CT_NEG1 > 0) / nrow(dfhypohyper),
                        percentHyper_G2_EG1 = sum(dfhypohyper$DiffMeth_CT_EG1 > 0) / nrow(dfhypohyper))
    
    return(result)
  }
  
  ## Loop over that for all BP:
  dfperchyper = data.frame(t(sapply(1:length(vecBP), getDFpercentHyper)))
  dfperchyper$BP = unlist(dfperchyper$BP)
  dfperchyper$percentHyper_G1=unlist(dfperchyper$percentHyper_G1)
  dfperchyper$percentHyper_G2_NEG1=unlist(dfperchyper$percentHyper_G2_NEG1)
  dfperchyper$percentHyper_G2_EG1=unlist(dfperchyper$percentHyper_G2_EG1)
  
  ggplot(dfperchyper)+
    geom_abline(slope = 1) +
    geom_label(aes(x=percentHyper_G2_NEG1, y=percentHyper_G1,label=BP), col = "black")+
    geom_label(aes(x=percentHyper_G2_EG1, y=percentHyper_G1,label=BP), col = "red")+
    xlab("Percentage of parDMS in G2 where methylation Exposed > methylation Control")+
    ylab("Percentage of parDMS in G1 where methylation Exposed > methylation Control")+
    xlim(0.25,0.75)+ylim(0.25,0.75)
}

# ## Manhattan plot of the 836!


## To do: separate into "get Differential methylation" and "analyse differential methylations"
