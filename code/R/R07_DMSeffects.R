# Each script sources the previous script of the pipeline if needed
source("R06_GlobalMethylationProfile.R")

## Load genome annotation
## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
## The default option is to take -1000,+1000bp around the TSS and you can change that. 
## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream
annotBed12=readTranscriptFeatures("../../gitignore/bigdata/06GynoAnnot/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",
                                  remove.unusual = FALSE, up.flank = 1500, down.flank = 500)
## Change recursively the gene names to keep only ID
getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
for (i in 1:length(annotBed12)){
  annotBed12[[i]]$name <- getName(annotBed12[[i]]$name)
}

# Load the calculated DMS by BP
DMBPlist <- readRDS("../../data/DiffMeth/DMperBP_list_dec2022.RDS")
DMSBPlist <- lapply(DMBPlist, "[[", 1)
DMRBPlist <- lapply(DMBPlist, "[[", 2)

## Description of DMS in the four different comparisons
# Differential methylation by brother pair (sex as covariate):
# 1.  CC-TC = CONTROL fish (parent CvsT)
# 2.  CT-TT = TREATMENT fish (parent CvsT)
# 3.  CC-CT = fish from CONTROL parents (G2 CvsT)
# 4.  TC-TT = fish from TREATMENT parents (G2 CvsT)

### Attribute plot of the DMS in every BP for the 4 comparisons
## All DMS are stored in DMSBPlist by brother pair

## Vector of all 4 comparisons
vecCompa <- c("CC_TC", "CT_TT", "CC_CT", "TC_TT")
vecCompaVerbose <- c("Control offspring in control vs infected parents", "Infected offspring in control vs infected parents", "Control vs infected offspring from control parent", "Control vs infected offspring from infected parent") # this is useful in plots

## Extract DMS for all 4 comparisons (by position)
myPosList = lapply(DMSBPlist, lapply, function(x){paste(x$chr, x$end)})

## Subselect those DMS present in at least 4 out of 8 BP
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

## Prepare df for complexUpset
getUpsetDF = function(i, NBP){ # for a given comparison
  A = get2keep(vecCompa[i], NBP)
  A2 = lapply(A, function(x){
    x = data.frame(x)    # vector of DMS as df
    names(x) = "DMS"    # name each CpG
    return(x)
  })
  ## Add BP name
  for (i in 1:length(names(A2))){
    A2[[i]]["BP"] = names(A2)[i]
  }
  # make a dataframe
  A2 = A2 %>% purrr::reduce(full_join, by = "DMS")
  # names column with BP id
  for (i in 2:ncol(A2)) {names(A2)[i] = unique(A2[!is.na(A2[i]), i])}
  # replace by 0 or 1 the DMS absence/presence
  A = data.frame(apply(A2[2:9], 2, function(x) ifelse(is.na(x), 0, 1)))
  # add DMS
  A$DMS = A2$DMS
  return(A)
}

#Make upset plots for DMS found in at least 6 BP:
for (i in 1:4){
  df = getUpsetDF(i, NBP = 6)
  print(ComplexUpset::upset(
    df,
    names(df)[1:8],
    width_ratio=0.1,
    themes=upset_default_themes(text=element_text(size=15)),
    sort_intersections_by=c('degree', 'cardinality'),
    queries=query_by_degree(
      df,  names(df)[1:8],
      params_by_degree=list(
        '1'=list(color='red', fill='red'),
        '2'=list(color='purple', fill='purple'),
        '3'=list(color='blue', fill='blue'),
        '4'=list(color='grey', fill='grey'),
        '5'=list(color='red', fill='red'),
        '6'=list(color='purple', fill='purple'),
        '7'=list(color='blue', fill='blue'),
        '8'=list(color='green', fill='green')
      ),
      only_components=c("intersections_matrix", "Intersection size")
    )) + ggtitle(paste0("Differentially methylated sites found in more than six brother pairs in the comparison: \n", vecCompaVerbose[i]))) #+ theme(plot.title = element_text(size = 30)))
}

## Statistical setup

# PARENTAL effect: DMS found in either CC-TC or CT-TT comparisons
# OFFSPRING effect: DMS found in either CC-CT or TC-TT comparisons
# INTERACTION effects: DMS found in CC-CT which show a differential methylation (not necessarily significant)
# in the opposite direction in TC-TT, or inversely

### Identify the different DMS groups
## Creates a data frame with BP, DMS, effect and comparison, when the DMS is found in 4 BP or more
get2keepFULLdfBP = function(Compa, NBP = 4, myeffect = "Paternal"){
  x <- lapply(myPosList, function(x){unlist(x[[paste0("DMS_15pc_BP_", Compa)]])})
  mylist = list() # empty list
  for (i in 1:length(names(x))){
    a=data.frame(BP=names(x)[i], DMS=x[[names(x)[i]]])
    mylist[[i]] = a
    names(mylist)[i] = names(x)[i]
  } # fill full list with df containing BP and DMS
  mydf = Reduce(function(...) merge(..., all=T), mylist) # concatenate in a df
  # Add information
  mydf$globalEffect = myeffect
  mydf$comparison = Compa
  # Keep the DMS present in 4 families minimum:
  tokeep = names(table(mydf$DMS)[table(mydf$DMS)>=NBP])
  mydf = mydf[mydf$DMS %in% tokeep,]
  return(mydf)
}

# Identify DMS in N BP or more (if needed to justify):
## Creates a data frame with BP, DMS, effect and comparison, when the DMS is found in N BP or more
## run over 1 to 8 BP

getDMS_runOver1to8BP <- function(NBP){
  # PARENTAL effect: DMS found in either CC-TC or CT-TT comparisons
  df_PaternalEffect = unique(rbind(get2keepFULLdfBP(Compa = vecCompa[1], NBP = NBP, myeffect = "G1all"),
                                   get2keepFULLdfBP(Compa = vecCompa[2], NBP = NBP, myeffect = "G1all")))
  DMS_PaternalEffect = unique(df_PaternalEffect$DMS)
  
  # OFFSPRING effect: DMS found in either CC-CT or TC-TT comparisons
  df_OffspringEffect = unique(rbind(get2keepFULLdfBP(Compa = vecCompa[3], NBP = NBP, myeffect = "G2all"),
                                    get2keepFULLdfBP(Compa = vecCompa[4], NBP = NBP, myeffect = "G2all")))
  DMS_OffspringEffect = unique(df_OffspringEffect$DMS)
  
  ##########
  ## OVERLAP=DMS found in G1 & G2
  # case 1 -> INTERACTION effects DMS found in CC-CT which show a differential methylation in the opposite direction in TC-TT, or inversely (reaction norm are inversed) in most of the brother pairs (5 or more)
  # case 2 -> ADDITIVE effect: no slope inversion (4 or more)
  
  # Get the raw methylation from the DMS which are in both Paternal and Offspring effects
  subUnite = methylKit::select(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                               which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$end) %in%
                                       intersect(DMS_OffspringEffect, DMS_PaternalEffect)))
  
  # Get mean methylation per brother pair, per treatment:
  getMeanMeth <- function(subUnite, BP, mytrt){
    metadata = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP & fullMetadata_OFFS$trtG1G2 %in% mytrt, ]
    myuniteCov = reorganize(methylObj = subUnite, treatment = metadata$trtG1G2_NUM, sample.ids = metadata$ID)
    ## remove bases where NO fish in this BP has a coverage
    myuniteCov = methylKit::select(myuniteCov, which(!is.na(rowSums(percMethylation(myuniteCov)))))
    # calculate mean methylation
    df = data.frame(DMS = paste(myuniteCov$chr, myuniteCov$end), meanMeth = rowMeans(percMethylation(myuniteCov)), trt = mytrt, BP = BP)
    return(df)
  }
  
  # We will apply the following function to all BP and all trt:
  vecBP <- unique(fullMetadata_OFFS$brotherPairID)
  vectrt <- unique(fullMetadata_OFFS$trtG1G2)
  
  ## Loop over all BP & trt
  df = data.frame(DMS=NULL, meanMeth=NULL, trt=NULL, BP=NULL) # empty df
  for (i in 1:length(vecBP)){
    for (j in 1:length(vectrt)){
      subdf = getMeanMeth(subUnite = subUnite, BP = vecBP[[i]], mytrt = vectrt[[j]])
      df = rbind(df, subdf)
    }
  }
  
  ## Add G1 & G2 trt
  df$G1trt = ifelse(df$trt %in% c("NE_control", "NE_exposed"), "control", "infected")
  df$G2trt = ifelse(df$trt %in% c("NE_control", "E_control"), "control", "infected")
  
  ## cut by G1 trt & merge
  dfcp = df[df$G1trt %in% "control"  ,]
  dfcpco = dfcp[dfcp$G2trt %in% "control",]; dfcpio = dfcp[dfcp$G2trt %in% "infected",]
  dfcp = merge(dfcpco, dfcpio, by = c("DMS", "BP")) %>%
    mutate(meanDiffMeth=meanMeth.y - meanMeth.x) %>% dplyr::select(c("DMS", "BP", "meanDiffMeth"))
  dfip = df[df$G1trt %in% "infected",]
  dfipco = dfip[dfip$G2trt %in% "control",]; dfipio = dfip[dfip$G2trt %in% "infected",]
  dfip = merge(dfipco, dfipio, by = c("DMS", "BP")) %>%
    mutate(meanDiffMeth=meanMeth.y - meanMeth.x) %>% dplyr::select(c("DMS", "BP", "meanDiffMeth"))
  df2=merge(dfcp,dfip, by=c("DMS", "BP"))
  
  ## interaction if slope inversion/additive if not
  df2$inversionSlopeReactionNorms = FALSE
  df2[!sign(df2$meanDiffMeth.x) == sign(df2$meanDiffMeth.y),"inversionSlopeReactionNorms"] = TRUE
  names(df2)[names(df2) %in% "meanDiffMeth.x"] = "meanDiffMeth.controlG1"
  names(df2)[names(df2) %in% "meanDiffMeth.y"] = "meanDiffMeth.infectedG1"
  
  ####################
  # Get a vector of DMS for each category:
  ## A DMS is "interaction" if there are more often slope inversion than not
  DMS_G1G2interactionEffect = df2 %>% group_by(DMS) %>% dplyr::summarise(count=n(), inversionSlopeRate=sum(inversionSlopeReactionNorms)/count,                                                                 
                                                                         effect=ifelse(inversionSlopeRate>0.5, "interaction", "additive")) %>%
    dplyr::filter(effect=="interaction") %>% dplyr::select(DMS) %>% unlist()
  
  DMS_G1G2additiveEffect = df2 %>% group_by(DMS) %>%
    dplyr::summarise(count=n(), inversionSlopeRate=sum(inversionSlopeReactionNorms)/count, effect=ifelse(inversionSlopeRate>0.5, "interaction", "additive")) %>%
    dplyr::filter(effect=="additive") %>% dplyr::select(DMS) %>% unlist()
  
  DMS_G1onlyEffect = DMS_PaternalEffect[!DMS_PaternalEffect %in% c(DMS_G1G2interactionEffect, DMS_G1G2additiveEffect)]
  
  DMS_G2onlyEffect = DMS_OffspringEffect[!DMS_OffspringEffect %in% c(DMS_G1G2interactionEffect, DMS_G1G2additiveEffect)]
  
  return(data.frame(NbrDMS = c(length(DMS_G1onlyEffect),
                               length(DMS_G2onlyEffect),
                               length(DMS_G1G2interactionEffect),
                               length(DMS_G1G2additiveEffect)), 
                    DMSgroup = c("Paternal only", "Offspring only", "interaction", "additive"),
                    NbrBP = rep(NBP, 4)))
}

A=do.call(rbind, lapply(1:6, getDMS_runOver1to8BP)) 

ggplot(A, aes(fill=DMSgroup, y=NbrDMS, x=NbrBP)) + 
  geom_bar(position="stack", stat="identity")+
  scale_y_log10()

A %>% head

# Identify DMS in 4 BP or more:

# PARENTAL effect: DMS found in either CC-TC or CT-TT comparisons
df_PaternalEffect_4BPmin = unique(rbind(get2keepFULLdfBP(Compa = vecCompa[1], NBP = 4, myeffect = "G1all"),
                                        get2keepFULLdfBP(Compa = vecCompa[2], NBP = 4, myeffect = "G1all")))
DMS_PaternalEffect_4BPmin = unique(df_PaternalEffect_4BPmin$DMS)

# OFFSPRING effect: DMS found in either CC-CT or TC-TT comparisons
df_OffspringEffect_4BPmin = unique(rbind(get2keepFULLdfBP(Compa = vecCompa[3], NBP = 4, myeffect = "G2all"),
                                         get2keepFULLdfBP(Compa = vecCompa[4], NBP = 4, myeffect = "G2all")))
DMS_OffspringEffect_4BPmin = unique(df_OffspringEffect_4BPmin$DMS)

##########
## OVERLAP=DMS found in G1 & G2
# case 1 -> INTERACTION effects DMS found in CC-CT which show a differential methylation in the opposite direction in TC-TT, or inversely (reaction norm are inversed) in most of the brother pairs (5 or more)
# case 2 -> ADDITIVE effect: no slope inversion (4 or more)

# Get the raw methylation from the DMS which are in both Paternal and Offspring effects
subUnite = methylKit::select(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                             which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$end) %in%
                                     intersect(DMS_OffspringEffect_4BPmin, DMS_PaternalEffect_4BPmin)))

# Get mean methylation per brother pair, per treatment:
getMeanMeth <- function(subUnite, BP, mytrt){
  metadata = fullMetadata_OFFS[fullMetadata_OFFS$brotherPairID %in% BP & fullMetadata_OFFS$trtG1G2 %in% mytrt, ]
  myuniteCov = reorganize(methylObj = subUnite, treatment = metadata$trtG1G2_NUM, sample.ids = metadata$ID)
  ## remove bases where NO fish in this BP has a coverage
  myuniteCov = methylKit::select(myuniteCov, which(!is.na(rowSums(percMethylation(myuniteCov)))))
  # calculate mean methylation
  df = data.frame(DMS = paste(myuniteCov$chr, myuniteCov$end), meanMeth = rowMeans(percMethylation(myuniteCov)), trt = mytrt, BP = BP)
  return(df)
}

# We will apply the following function to all BP and all trt:
vecBP <- unique(fullMetadata_OFFS$brotherPairID)
vectrt <- unique(fullMetadata_OFFS$trtG1G2)

## Loop over all BP & trt
df = data.frame(DMS=NULL, meanMeth=NULL, trt=NULL, BP=NULL) # empty df
for (i in 1:length(vecBP)){
  for (j in 1:length(vectrt)){
    subdf = getMeanMeth(subUnite = subUnite, BP = vecBP[[i]], mytrt = vectrt[[j]])
    df = rbind(df, subdf)
  }
}

## Add G1 & G2 trt
df$G1trt = ifelse(df$trt %in% c("NE_control", "NE_exposed"), "control", "infected")
df$G2trt = ifelse(df$trt %in% c("NE_control", "E_control"), "control", "infected")

## cut by G1 trt & merge
dfcp = df[df$G1trt %in% "control"  ,]
dfcpco = dfcp[dfcp$G2trt %in% "control",]; dfcpio = dfcp[dfcp$G2trt %in% "infected",]
dfcp = merge(dfcpco, dfcpio, by = c("DMS", "BP")) %>%
  mutate(meanDiffMeth=meanMeth.y - meanMeth.x) %>% dplyr::select(c("DMS", "BP", "meanDiffMeth"))
dfip = df[df$G1trt %in% "infected",]
dfipco = dfip[dfip$G2trt %in% "control",]; dfipio = dfip[dfip$G2trt %in% "infected",]
dfip = merge(dfipco, dfipio, by = c("DMS", "BP")) %>%
  mutate(meanDiffMeth=meanMeth.y - meanMeth.x) %>% dplyr::select(c("DMS", "BP", "meanDiffMeth"))
df2=merge(dfcp,dfip, by=c("DMS", "BP"))

## interaction if slope inversion/additive if not
df2$inversionSlopeReactionNorms = FALSE
df2[!sign(df2$meanDiffMeth.x) == sign(df2$meanDiffMeth.y),"inversionSlopeReactionNorms"] = TRUE
names(df2)[names(df2) %in% "meanDiffMeth.x"] = "meanDiffMeth.controlG1"
names(df2)[names(df2) %in% "meanDiffMeth.y"] = "meanDiffMeth.infectedG1"

####################
# Get a vector of DMS for each category:
## A DMS is "interaction" if there are more often slope inversion than not
DMS_G1G2interactionEffect_4BPmin = df2 %>% group_by(DMS) %>% dplyr::summarise(count=n(), inversionSlopeRate=sum(inversionSlopeReactionNorms)/count,                                                                 effect=ifelse(inversionSlopeRate>0.5, "interaction", "additive")) %>%
  dplyr::filter(effect=="interaction") %>% dplyr::select(DMS) %>% unlist()

DMS_G1G2additiveEffect_4BPmin = df2 %>% group_by(DMS) %>%
  dplyr::summarise(count=n(), inversionSlopeRate=sum(inversionSlopeReactionNorms)/count, effect=ifelse(inversionSlopeRate>0.5, "interaction", "additive")) %>%
  dplyr::filter(effect=="additive") %>% dplyr::select(DMS) %>% unlist()

DMS_G1onlyEffect_4BPmin = DMS_PaternalEffect_4BPmin[!DMS_PaternalEffect_4BPmin %in% c(DMS_G1G2interactionEffect_4BPmin, DMS_G1G2additiveEffect_4BPmin)]

DMS_G2onlyEffect_4BPmin = DMS_OffspringEffect_4BPmin[!DMS_OffspringEffect_4BPmin %in% c(DMS_G1G2interactionEffect_4BPmin, DMS_G1G2additiveEffect_4BPmin)]

####################
# Make a BIG df with all DMS, effects and BP (this time, the effects are DMS-BP specific, not global)
df_effects_full = merge(unique(df_PaternalEffect_4BPmin[c("BP", "DMS","globalEffect")]),
                        unique(df_OffspringEffect_4BPmin[c("BP", "DMS","globalEffect")]), by=c("BP", "DMS"), all=T)

df_effects_full = merge(df_effects_full, df2, all=T)
df_effects_full$effectBPlevel[df_effects_full$globalEffect.x == "G1all"] = "G1"
df_effects_full$effectBPlevel[df_effects_full$globalEffect.y == "G2all"] = "G2"
df_effects_full$effectBPlevel[df_effects_full$inversionSlopeReactionNorms == TRUE & 
                                df_effects_full$globalEffect.x == "G1all" & df_effects_full$globalEffect.y == "G2all"] = "inter"
df_effects_full$effectBPlevel[df_effects_full$inversionSlopeReactionNorms == FALSE & 
                                df_effects_full$globalEffect.x == "G1all" & df_effects_full$globalEffect.y == "G2all"]= "addit"

#rm junk
rm(subUnite, df, dfcp, dfcpco, dfcpio, dfip, dfipco, dfipio, df2)

# Plot a Venn diagram
ggVennDiagram(list("Paternal effect" = DMS_PaternalEffect_4BPmin, "Offspring effect" = DMS_OffspringEffect_4BPmin, "InteractionEffects" = DMS_G1G2interactionEffect_4BPmin),
              label_alpha = 0) + scale_color_manual(values = c(1,1,1))+
  scale_fill_gradient(low="white",high = "yellow") + theme(legend.position = "none")

# Save:
pdf(file = "../../dataOut/DMS3groupsVenn.pdf", width = 7, height = 6)
ggVennDiagram(list("Paternal effect" = DMS_PaternalEffect_4BPmin, "Offspring effect" = DMS_OffspringEffect_4BPmin, "InteractionEffects" = DMS_G1G2interactionEffect_4BPmin),
              label_alpha = 0) + scale_color_manual(values = c(1,1,1))+
  scale_fill_gradient(low="white",high = "yellow") + theme(legend.position = "none")
dev.off()

#######################
## Where are these DMS?
DMSvec=unique(c(DMS_G1onlyEffect_4BPmin, DMS_G2onlyEffect_4BPmin, DMS_G1G2additiveEffect_4BPmin, DMS_G1G2interactionEffect_4BPmin))

getFeature <- function(DMSvec){
  # Change the DMS vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=sapply(strsplit(DMSvec, " "), `[`, 1), 
                                                  start=sapply(strsplit(DMSvec, " "), `[`, 2),
                                                  end=sapply(strsplit(DMSvec, " "), `[`, 2),
                                                  DMS=DMSvec), keep.extra.columns = T)
  annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = annotBed12)
}

A=getFeature(DMSvec = DMSvec)
print(paste0("Positions of the ", length(DMSvec)," DMS:"))
print(A)

# A1 = getFeature(DMSvec = DMS_G1onlyEffect_4BPmin)
# print(paste0("Positions of the ", length(DMS_G1onlyEffect_4BPmin)," intergenerational DMS:"))
# print(A1@precedence)
# 
# A2 = getFeature(DMSvec = DMS_G2onlyEffect_4BPmin)
# print(paste0("Positions of the ", length(DMS_G1onlyEffect_4BPmin)," infection-induced DMS:"))
# print(A2@precedence)
# 
# A3 = getFeature(DMSvec = DMS_G1G2additiveEffect_4BPmin)
# print(paste0("Positions of the ", length(DMS_G1onlyEffect_4BPmin)," additive DMS:"))
# print(A3@precedence)
# 
# A4 = getFeature(DMSvec = DMS_G1G2interactionEffect_4BPmin)
# print(paste0("Positions of the ", length(DMS_G1onlyEffect_4BPmin)," interaction DMS:"))
# print(A4@precedence)

#######################
## Are the positions of DMS on features random? Comparison with sequenced CpGs which are not DMS
# A=getFeature(DMSvec = DMSvec)
AnonDMS= getFeature(DMSvec = paste(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$start))

ChiTable1 = merge((A@members %>% data.frame() %>% mutate(feature=ifelse(prom==1, "promoter", 
                                                                        ifelse(exon==1, "exon",
                                                                               ifelse(intron==1, "intron", "intergenic")))) %>% 
                     dplyr::select(feature) %>% 
                     table %>% melt %>% dplyr::rename(DMS=value)),
                  (AnonDMS@members %>% data.frame() %>% mutate(feature=ifelse(prom==1, "promoter", 
                                                                              ifelse(exon==1, "exon",
                                                                                     ifelse(intron==1, "intron", "intergenic")))) %>% 
                     dplyr::select(feature) %>% 
                     table %>% melt %>% dplyr::rename(nonDMS=value)))

chisq.test(ChiTable1$DMS, ChiTable1$nonDMS)
# X-squared = 12, df = 9, p-value = 0.2133

donutDF = ChiTable1 %>% dplyr::mutate(percDMS=DMS/sum(DMS)*100, percNonDMS=nonDMS/sum(nonDMS)*100) %>%
  dplyr::select("feature", "percDMS", "percNonDMS") %>% melt

ggplot(donutDF, aes(x = variable, y = value, fill = feature)) +
  geom_col() +  scale_x_discrete(limits = c(" ", "percNonDMS","percDMS")) +
  scale_fill_viridis_d()+
  coord_polar("y")

## Positions on chromosomes?
df = A@members %>% data.frame() %>% mutate(feature=ifelse(prom==1, "promoter", 
                                                          ifelse(exon==1, "exon",
                                                                 ifelse(intron==1, "intron", "intergenic"))))%>%
  
  mutate(DMS=DMSvec, 
         chr=sapply(strsplit(DMSvec, " "), `[`, 1)%>% str_remove(.,"Gy_chr"))

melt(table(df$feature))

df2 = melt(table(df$feature, df$chr))
names(df2)=c("feature","chromosome", "nDMS")

# Number of positions sequenced on each chromosome
df3 = table(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr) %>% data.frame() %>% 
  mutate(chromosome=str_remove(Var1,"Gy_chr"), nCpG=Freq) 

df=merge(df2, df3) %>% dplyr::select(c(feature, chromosome, nDMS, nCpG))%>% 
  mutate(percent=nDMS/nCpG)

df$feature=factor(df$feature, levels = c("promoter", "exon", "intron", "intergenic"))

pdf(file = "../../dataOut/suppl1_featureDist.pdf", width = 8, height = 6)
ggplot(df, aes(x=chromosome, y=percent, fill=feature))+
  geom_bar(stat = "identity")+
  ylab("Percentage of DMS among CpG sequenced")+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values = c("red", "purple", "blue", "grey"))
dev.off()

dfChi = df %>% group_by(chromosome) %>% dplyr::summarise(nDMS=sum(nDMS)) %>%
  merge(df3) 

chisq.test(dfChi[c("nDMS", "nCpG")])
# X-squared = 233.87, df = 19, p-value < 2.2e-16

dfChi %>% mutate(percent=nDMS/(nCpG)*100) %>% arrange(percent) 
# range from 0.13% (XV) to 0.39% of CpG beign DMS (XVIII)
