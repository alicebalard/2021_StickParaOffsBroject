getInteractionDMS <- function(){
  
  ####### Part 1: first direction
  ## Extract a data frame with DMS and meth.diff for one given comparison (vecCompa[3]= CC_CT):
  mydms_CC_CT = lapply(
    lapply(DMSBPlist, lapply, function(x){data.frame(DMS = paste(x$chr, x$end), meth.diff = x$meth.diff)}), 
    function(x){x[[paste0("DMS_15pc_BP_", vecCompa[3])]]})
  
  ## Add BP in the name of the column containing meth.diff:
  for (i in 1:length(names(mydms_CC_CT))){
    names(mydms_CC_CT[[i]])[names(mydms_CC_CT[[i]]) %in% "meth.diff"] = paste0("meth.diff_", names(mydms_CC_CT)[i])
  }
  
  # merge all elements of the list in one big data frame
  mydms_CC_CT = Reduce(function(...) merge(..., all=T), mydms_CC_CT)
  
  # keep DMS found in at least 4 BP:
  mydms_CC_CT_4BPmin = mydms_CC_CT[rowSums(is.na(mydms_CC_CT[2:ncol(mydms_CC_CT)]))< 4, ]
  
  # Subselect the original unite object for these DMS only 
  subUnite = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, 
                               which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, uniteCov14_G2_woSexAndUnknowChrOVERLAP$end) %in% mydms_CC_CT_4BPmin$DMS))
  
  ### Per brother pair:
  smallgetdiffmeth <- function(BP, mytrt = c("E_control", "E_exposed")){# by default in (TC_TT)
    # Subselect for the samples of the other offspring effect comparison & a given BP
    metadata = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP & fullMetadata_OFFS$trtG1G2 %in% mytrt, ]
    myuniteCov = reorganize(methylObj = subUnite, treatment = metadata$trtG1G2_NUM, sample.ids = metadata$ID)
    ## remove bases where NO fish in this BP has a coverage
    myuniteCov = methylKit::select(myuniteCov, which(!is.na(rowSums(percMethylation(myuniteCov)))))
    # calculate differential methylation (without statistical significance) in each brother pair
    if (length(table(metadata$Sex)) == 1){ 
      myDiffMeth=calculateDiffMeth(myuniteCov)
    } else { # if more than 1 sex, we add a covariate
      cov = data.frame(Sex = metadata$Sex)
      myDiffMeth=calculateDiffMeth(myuniteCov, covariates = cov)
    } 
    return(myDiffMeth)
  }
  
  # We will apply the following function to all BP:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  
  ## Loop over all BP
  DMlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMlist[[i]] <- smallgetdiffmeth(BP = vecBP[[i]])
  } 
  names(DMlist) <- vecBP
  
  DMlist = lapply(DMlist, function(x){data.frame(DMS = paste(x$chr, x$end), meth.diff = x$meth.diff)})
  
  ## Add BP in the name of the column containing meth.diff:
  for (i in 1:length(names(DMlist))){
    names(DMlist[[i]])[names(DMlist[[i]]) %in% "meth.diff"] = paste0("meth.diff.otherComp_", names(DMlist)[i])
  }
  
  # merge all elements of the list in one big data frame
  DMlist = Reduce(function(...) merge(..., all=T), DMlist)
  
  fullDF = data.frame(DMS=NULL)
  for (BP in vecBP){
    df=merge(DMlist[c("DMS", paste0("meth.diff.otherComp_", BP))], mydms_CC_CT_4BPmin[c("DMS", paste0("meth.diff_", BP))])
    df$BP = BP
    names(df) = gsub(paste0("_", BP), "", names(df))
    melt(df)
    # Merge with the rest of BP df
    fullDF=merge(fullDF, df, all = T)
  }
  fullDF = na.omit(fullDF)
  
  # compare the sign of the slope in both cases:
  fullDF$sign.meth.diff = sign(fullDF$meth.diff)
  fullDF$sign.meth.diff.otherComp = sign(fullDF$meth.diff.otherComp)
  
  fullDF$compSlope = fullDF$sign.meth.diff != fullDF$sign.meth.diff.otherComp
  
  # We calculate mean(compSlope) = the percentage of BP in which the slope is in the same direction, for each DMS
  fullDFSum = fullDF %>% group_by(DMS) %>% dplyr::summarise(mean = mean(compSlope), n=n()) %>% data.frame()
  
  # We select the DMS for which in more than 50% of the BP considered, the slope is INVERSED in both offspring effect comparisons
  set1DMSInteractions = fullDFSum$DMS[fullDFSum$mean >0.5]
  
  #####################
  ####### Part 2: second direction
  ## Extract a data frame with DMS and meth.diff for one given comparison (vecCompa[4]= TC_TT):
  mydms_TC_TT = lapply(
    lapply(DMSBPlist, lapply, function(x){data.frame(DMS = paste(x$chr, x$end), meth.diff = x$meth.diff)}), 
    function(x){x[[paste0("DMS_15pc_BP_", vecCompa[4])]]})
  
  ## Add BP in the name of the column containing meth.diff:
  for (i in 1:length(names(mydms_TC_TT))){
    names(mydms_TC_TT[[i]])[names(mydms_TC_TT[[i]]) %in% "meth.diff"] = paste0("meth.diff_", names(mydms_TC_TT)[i])
  }
  
  # merge all elements of the list in one big data frame
  mydms_TC_TT = Reduce(function(...) merge(..., all=T), mydms_TC_TT)
  
  # keep DMS found in at least 4 BP:
  mydms_TC_TT_4BPmin = mydms_TC_TT[rowSums(is.na(mydms_TC_TT[2:ncol(mydms_TC_TT)]))< 4, ]
  
  # Subselect the original unite object for these DMS only 
  subUnite = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, 
                               which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, uniteCov14_G2_woSexAndUnknowChrOVERLAP$end) %in% mydms_TC_TT_4BPmin$DMS))
  
  ### Per brother pair:
  # We will apply the following function to all BP:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  
  ## Loop over all BP
  DMlist <- list() # empty plot list
  for (i in 1:length(vecBP)){
    DMlist[[i]] <- smallgetdiffmeth(BP = vecBP[[i]], mytrt = c("NE_control", "NE_exposed"))
  } 
  names(DMlist) <- vecBP
  
  DMlist = lapply(DMlist, function(x){data.frame(DMS = paste(x$chr, x$end), meth.diff = x$meth.diff)})
  
  ## Add BP in the name of the column containing meth.diff:
  for (i in 1:length(names(DMlist))){
    names(DMlist[[i]])[names(DMlist[[i]]) %in% "meth.diff"] = paste0("meth.diff.otherComp_", names(DMlist)[i])
  }
  
  # merge all elements of the list in one big data frame
  DMlist = Reduce(function(...) merge(..., all=T), DMlist)
  
  fullDF = data.frame(DMS=NULL)
  for (BP in vecBP){
    df=merge(DMlist[c("DMS", paste0("meth.diff.otherComp_", BP))], mydms_TC_TT_4BPmin[c("DMS", paste0("meth.diff_", BP))])
    df$BP = BP
    names(df) = gsub(paste0("_", BP), "", names(df))
    melt(df)
    # Merge with the rest of BP df
    fullDF=merge(fullDF, df, all = T)
  }
  fullDF = na.omit(fullDF)
  
  # compare the sign of the slope in both cases:
  fullDF$sign.meth.diff = sign(fullDF$meth.diff)
  fullDF$sign.meth.diff.otherComp = sign(fullDF$meth.diff.otherComp)
  
  fullDF$compSlope = fullDF$sign.meth.diff != fullDF$sign.meth.diff.otherComp
  
  # We calculate mean(compSlope) = the percentage of BP in which the slope is in the same direction, for each DMS
  fullDFSum = fullDF %>% group_by(DMS) %>% dplyr::summarise(mean = mean(compSlope), n=n()) %>% data.frame()
  
  # We select the DMS for which in more than 50% of the BP considered, the slope is INVERSED in both offspring effect comparisons
  set2DMSInteractions = fullDFSum$DMS[fullDFSum$mean >0.5]
  
  setDMSInteractions = unique(c(set1DMSInteractions, set2DMSInteractions))
  
  return(setDMSInteractions)
}

setDMSInteractions = getInteractionDMS()
