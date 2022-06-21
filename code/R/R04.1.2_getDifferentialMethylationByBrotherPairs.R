## Differential methylation analyses BY BROTHER PAIR
## A. Balard
## February 2022

machine="apocrita" # define the machine we work on
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

############################
## II. Bottom-Up approach ##
############################

#######################################
## Attribute plot of DMS in every BP ##
#######################################

nrow(uniteCov14_G2_woSexAndUnknowChrOVERLAP) # methylBase object with 1001880 rows

## Calculate DMS between G2 from G1C & G2 from G1T, by brother pair
# getDiffMeth(uniteCov14_G2_woSexAndUnknowChrOVERLAP, fullMetadata_PAR)

uniteCov14_G2_woSexAndUnknowChrOVERLAP

table(fullMetadata_OFFS$trtG1G2, fullMetadata_OFFS$trtG1G2_NUM)

###################################################################################################
# Calculate DMS CC-TC and TC-TT in all brother pairs (same treatment, different PARENTAL treatment)

getDMSperBP <- function(BP){
  ## Unite object for one Brother Pair:
  metadataBP_control = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                           fullMetadata_OFFS$trtG1G2 %in% c("NE_control", "E_control"), ]
  metadataBP_treatment = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP &
                                             fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"), ]
  
  ## Make 2 separate uniteCov, for each offspring trt:
  myuniteCovBP_control = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                    treatment = metadataBP_control$trtG1G2_NUM, sample.ids = metadataBP_control$ID)
  myuniteCovBP_treatment = reorganize(methylObj = uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                                      treatment = metadataBP_treatment$trtG1G2_NUM, sample.ids = metadataBP_treatment$ID)
  
  # remove bases where NO fish in this BP has a coverage
  myuniteCovBP_control = methylKit::select(myuniteCovBP_control, which(!is.na(rowSums(percMethylation(myuniteCovBP_control)))))
  myuniteCovBP_treatment = methylKit::select(myuniteCovBP_treatment, which(!is.na(rowSums(percMethylation(myuniteCovBP_treatment)))))
  
  # Calculate differential methylation
  myDiffMethBP = calculateDiffMeth(myuniteCovBP, mc.cores = 10)
  
  myDMS_15pc_BP = getMethylDiff(myDiffMethBP, difference=15, qvalue=0.01)
  
  # We select the bases that have q-value<0.01 and percent methylation difference larger than 15%, sex as covariate
  DMS_15pc_BP_control = getDiffMeth(myuniteCov = myuniteCovBP_control, myMetadata = metadataBP_control, mccores = 10, mydif = 15)
  DMS_15pc_BP_treatment = getDiffMeth(myuniteCov = myuniteCovBP_treatment, myMetadata = metadataBP_treatment, mccores = 10, mydif = 15)
  
  return(list(DMS_15pc_BP_control = DMS_15pc_BP_control, DMS_15pc_BP_treatment = DMS_15pc_BP_treatment))
}

# myPos = paste(myDMS_15pc_BP$chr, myDMS_15pc_BP$end)


run=TRUE
if (run==TRUE){
  ## Loop over all BP
  DMSlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMSlist[[i]] <- getDMSperBP(BP = vecBP[[i]])
  } 
  names(DMSlist) <- vecBP
  saveRDS(DMSlist, "./Rdata/DMS_BP_G2_list.RDS")
}

# DMS_BP_G2_CparvsTpar <- readRDS("Rdata/DMS_BP_G2_CparvsTpar.RDS")
# 
# ## Add fam to bp
# names(DMS_BP_G2_CparvsTpar) <- paste0(join(data.frame(BP = names(DMS_BP_G2_CparvsTpar)),
#                                            unique(data.frame(BP = fullMetadata_OFFS$brotherPairID, fam = fullMetadata_OFFS$Family)))[[2]],
#                                       "_", names(DMS_BP_G2_CparvsTpar))
# 
# library(UpSetR)
# 
# ## Find DMS present in at least 4 BP (half):
# x <- unlist(DMS_BP_G2_CparvsTpar)
# f <- table(x)
# tokeep <- names(f)[f >= 4] 
# DMS_BP_G2_CparvsTpar_INTER4 <- lapply(DMS_BP_G2_CparvsTpar, function(x){x[x %in% tokeep]})
# 
# ## Reorder by family:
# DMS_BP_G2_CparvsTpar_INTER4 <- DMS_BP_G2_CparvsTpar_INTER4[
#   names(DMS_BP_G2_CparvsTpar_INTER4)[order(names(DMS_BP_G2_CparvsTpar_INTER4))]]
# 
# # The majority of differences are within BP or family 
# UpSetR::upset(fromList(DMS_BP_G2_CparvsTpar_INTER4), order.by = "freq", nsets = 8, keep.order = T, 
#               sets = names(DMS_BP_G2_CparvsTpar_INTER4), shade.color = "grey")
# 
# ###################
# # Get annotation of the DMS present in at least 4 BP
# DMSvec <- unique(unlist(DMS_BP_G2_CparvsTpar_INTER4)) # N = 836
# 
# # Annotate that:
# # Change the vector into a methobject:
# df <- data.frame(chr=sapply(strsplit(DMSvec, " "), `[`, 1),
#                  start=sapply(strsplit(DMSvec, " "), `[`, 2), 
#                  end=sapply(strsplit(DMSvec, " "), `[`, 2))
# # get annotation
# anot <- getAnnotationFun(makeGRangesFromDataFrame(df))
# #356 genes
# 
# ggplot(anot, aes(x=start, y = nCpG)) + geom_point(alpha=.4) +
#   facet_grid(.~seqnames) + 
#   # geom_label(data = anot[anot$nCpG > 10,], aes(label = Note)) +
#   theme(axis.text.x=element_blank()) +
#   xlab("Genes with DMS present in at least 4 brother pairs")
# 
# ## let's find out about the top 3, on chr3, chr14 and chr21
# anot[anot$nCpG >11,]
# # Gy_chrXIV   27332   34913  7582 gasAcul20078-RA 
# # Similar to Ptpn11: Tyrosine-protein phosphatase non-receptor type 11 (Rattus norvegicus OX=10116)
# 
# # Gy_chrXXI 2580241 2584755  4515 gasAcul22312-RA 
# # Similar to Aoc1: Amiloride-sensitive amine oxidase [copper-containing] (Rattus norvegicus OX=10116)
# # Catalyzes the degradation of compounds such as putrescine, histamine, spermine, and spermidine, 
# # substances involved in allergic and immune responses, cell proliferation, tissue differentiation, tumor formation, and possibly apoptosis.
# 
# # Gy_chrIII 1700100 1716803 16704 gasAcul16782-RA 
# # Similar to KIF1A: Kinesin-like protein KIF1A (Homo sapiens OX=9606)
# # Motor for anterograde axonal transport of synaptic vesicle precursors. Also required for neuronal 
# # dense core vesicles (DCVs) transport to the dendritic spines and axons. The interaction calcium-dependent 
# # with CALM1 increases vesicle motility and interaction with the scaffolding proteins PPFIA2 and TANC2 recruits DCVs to synaptic sites.
# 
# ###################
# ## Some GO:
# 
# ## Gene universe: all genes which are present in the dataset
# 
# ## Gene sub-universe: all genes in DMS
# # in the GOterm universe we only want genes for which (1) CpGs were covered, and (2) which have a GOterm attributed:
# gene_universe <- data.frame(
#   subsetByOverlaps(GRanges(annotGff3), GRanges(uniteCov14_G2_woSexAndUnknowChrOVERLAP))) %>% # subselect covered CpGs
#   filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
#   filter(type %in% "gene")  %>% # keep all the 7416 genes with GO terms
#   dplyr::select(c("Name", "Ontology_term")) %>% 
#   mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, it's from Interproscan after Maker. Sticklebacks (biomart) have 82701 IEA and 63 ISS.
#   relocate("Ontology_term","go_linkage_type","Name") %>% 
#   unnest(Ontology_term) %>% # one GO per line (was a list before in this column)
#   data.frame()
# 
# # Create gene set collection
# goFrame <- GOFrame(gene_universe, organism="Gasterosteus aculeatus")
# goAllFrame <- GOAllFrame(goFrame)
# gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())
# 
# ## Create subuniverse: DMS between G2 from G1C & G2 from G1T, by brother pair, present in AT LEAST 4 BP/8
# sub_universe <- gene_universe %>%
#   subset(gene_universe$Name %in% unlist(anot$Parent))
# 
# ###### IMPORTANT NOTE from Mel: why conditional hypergeometric test?
# #The GO ontology is set up as a directed acyclic graph, where a parent term is comprised of all its child terms. If you do a standard
# #hypergeometric, you might e.g., find 'positive regulation of kinase activity' to be significant.
# #If you then test 'positive regulation of catalytic activity', which is a parent term, then it might be significant as well, but only because of
# #the terms coming from positive regulation of kinase activity.
# 
# #The conditional hypergeometric takes this into account, and only uses
# #those terms that were not already significant when testing a higher
# #order (parent) term.
# #####################
# 
# ## Run conditional hypergeometric test:
# runTestHypGeom <- function(sub_universe, onto){
#   ## Constructing a GOHyperGParams objects or KEGGHyperGParams objects from a GeneSetCollection
#   params_sequenced <- GSEAGOHyperGParams(name="GO_set",
#                                          geneSetCollection = gsc_universe,
#                                          geneIds = as.vector(unique(sub_universe[["Name"]])), # gene ids for the selected gene set
#                                          universeGeneIds = unique(gene_universe$Name),
#                                          ontology = onto, # A string specifying the GO ontology to use. Must be one of "BP", "CC", or "MF". (used with GO only)
#                                          pvalueCutoff = 0.05,
#                                          conditional = TRUE, # see note above
#                                          testDirection = "over") # over represented GO terms
#   
#   ## Run hypergeometric test:
#   return(hyperGTest(params_sequenced))
# }
# 
# ## Molecular functions
# GO_MF <- runTestHypGeom(sub_universe = sub_universe, onto = "MF")
# # Gene to GO MF Conditional test for over-representation 
# # 188 GO MF ids tested (17 have p < 0.05)
# # Selected gene set size: 151 
# # Gene universe size: 6397 
# GO_MF %>% summary() %>% head()
# 
# GO_CC <- runTestHypGeom(sub_universe = sub_universe, onto = "CC")
# # 67 GO CC ids tested (1 have p < 0.05)
# # Selected gene set size: 51 
# # Gene universe size: 2077 
# 
# GO_BP <- runTestHypGeom(sub_universe = sub_universe, onto = "BP")
# # 290 GO BP ids tested (10 have p < 0.05)
# # Selected gene set size: 98 
# # Gene universe size: 3577 
# 
# ## To save:
# # htmlReport(GO_MF_0633, file="../Routput/GO/GO_MF_0633.html")
# # write.table(as.data.frame(summary(GO_MF_0633)),"../Routput/GO/GO_MF_0633.txt", sep="\t", row.names = F)
# 
# ## Merge the df MP and BP
# A = GO_MF %>% summary %>% mutate(GOID = GOMFID, Type = "MF")
# B = GO_BP %>% summary %>% mutate(GOID = GOBPID, Type = "BP")
# C = GO_CC %>% summary %>% mutate(GOID = GOCCID, Type = "CC")
# 
# dfGO = rbind(A[-1], B[-1], C[-1])
# 
# plotGO = dfGO %>% ggplot(aes(x=Count/Size, y = Term)) +
#   geom_point(aes(color = Pvalue, size = Count)) +
#   scale_color_gradient(name="P value", low = "red", high = "blue") +
#   scale_size_continuous(name = "Gene number", range = c(1, 10), breaks = c(1,5,15,25))+
#   theme_bw() +
#   facet_grid(Type~.,scales="free",space = "free")
# 
# # pdf(xxx, width = 15, height = 20)
# plotGO
# # dev.off()
# dfGO
# 
# 
# ## Manhattan plot of the 836!
