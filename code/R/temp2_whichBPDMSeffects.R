##########
## OVERLAP=DMS found as G1 & G2 in at least 4 BP
# INTERACTION effects DMS found in CC-CT which show a differential methylation in the opposite direction in TC-TT, or inversely (reaction norm are inversed) 
# ADDITIVE effect: no slope inversion

# Get the raw methylation from the DMS which are in both Paternal and Offspring effects
subUnite = methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP, 
                             which(paste(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr, uniteCov14_G2_woSexAndUnknowChrOVERLAP$end) %in% 
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

# Keep DMS for which there is G1 AND G2 effect in at least 4 BP
df2=df2[paste(df2$BP, df2$DMS) %in% paste(df_PaternalEffect_4BPmin$BP, df_PaternalEffect_4BPmin$DMS) &
          paste(df2$BP, df2$DMS) %in% paste(df_OffspringEffect_4BPmin$BP, df_OffspringEffect_4BPmin$DMS),]
df2=df2[df2$DMS %in% names(table(df2$DMS)[table(df2$DMS)>=4]),]

# Keep rows for which there is an inversion of sign as INTERACTION, and no inversion as ADDITIVE
dfINTER=df2[!sign(df2$meanDiffMeth.x) == sign(df2$meanDiffMeth.y),] %>% dplyr::select(DMS, BP) %>% mutate(effect="G1:G2")
dfADDIT=df2[sign(df2$meanDiffMeth.x) == sign(df2$meanDiffMeth.y),] %>% dplyr::select(DMS, BP) %>% mutate(effect="G1+G2")

####################
# Make a BIG df with all DMS, effects and BP
df_effects_full = merge(unique(df_PaternalEffect_4BPmin[c("BP", "DMS","globalEffect")]),
                        unique(df_OffspringEffect_4BPmin[c("BP", "DMS","globalEffect")]), by=c("BP", "DMS"), all=T)

## Add interaction or additive (for each, in at least 4 BP)
df_effects_full = merge(df_effects_full, merge(dfINTER, dfADDIT, all = T), all = T)

## Add only G1 or only G2
df_effects_full$effect[df_effects_full$globalEffect.x %in% "G1all" & !df_effects_full$globalEffect.y %in% "G2all"] = "G1"
df_effects_full$effect[!df_effects_full$globalEffect.x %in% "G1all" & df_effects_full$globalEffect.y %in% "G2all"] = "G2"

## Add G1 or G2 (the overlap is only in some brother pairs)
df_effects_full$effect[is.na(df_effects_full$effect)] = "G1orG2"

#rm junk
rm(subUnite, df, dfcp, dfcpco, dfcpio, dfip, dfipco, dfipio, df2)

# check
hist(table(df_effects_full$DMS[df_effects_full$effect %in% c("G1+G2", "G1:G2")]))
hist(table(df_effects_full$DMS[df_effects_full$effect %in% c("G1orG2")]))

####################
# Get a vector of DMS for each category: WRONG

# Additive DMS: more often additive than  Interaction (and vice versa):
a=df_effects_full[df_effects_full$effect %in% "G1:G2",c("BP", "DMS")]
b=df_effects_full[df_effects_full$effect %in% "G1+G2",c("BP", "DMS")]

# Is considered "interaction" a DMS for which there is more often an inversion of slope than not
c = unique(c(a$DMS, b$DMS))
DMS_G1G2interactionEffect_4BPmin <- vector()
for (i in c){
  if(nrow(a[a$DMS %in% i,]) > nrow(b[b$DMS %in% i,])){
    DMS_G1G2interactionEffect_4BPmin <- c(intervec, i)
  }
}
DMS_G1G2additiveEffect_4BPmin = c[!c %in% DMS_G1G2interactionEffect_4BPmin]

## Assign DMS found ONLY in G1 or ONLY in G2 
DMS_G1onlyEffect_4BPmin = unique(df_effects_full$DMS[df_effects_full$effect %in% "G1"])[
  !unique(df_effects_full$DMS[df_effects_full$effect %in% "G1"]) %in% c]

DMS_G2onlyEffect_4BPmin = unique(df_effects_full$DMS[df_effects_full$effect %in% "G2"])[
  !unique(df_effects_full$DMS[df_effects_full$effect %in% "G2"]) %in% c]


# DMS_G1orG2Effect_4BPmin = unique(df_effects_full$DMS[df_effects_full$effect %in% "G1orG2"])

ggVennDiagram(list("G1"=DMS_G1onlyEffect_4BPmin, "G2"=DMS_G2onlyEffect_4BPmin, #"G1orG2"=DMS_G1orG2Effect_4BPmin, 
                   "G1+G2"=DMS_G1G2additiveEffect_4BPmin, "G1:G2"=DMS_G1G2interactionEffect_4BPmin),
              label_alpha = 0) + scale_fill_gradient(low="white",high = "yellow") 

DMS_G1onlyEffect_4BPmin

```

```{r VENNEffectDMS}
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