# Each script sources the previous script of the pipeline if needed
source("R10_intergenerationalEffect.R")

message("R11 starting...\n")

###################################################
## Correlation between methylation and phenotype ##
###################################################

## Correlation between methylation (after PCA) and phenotype (nbr worms, BCI)

# Approach:
# 1. Extract methylation values: raw beta values at DMS
# 2. PCA
# 3. Extract axes 1 & 2
# 4. Correlation with parasite load/BCI
# 5. If result positive, annotate the CpG associated with significant axis

### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values
DMS_allDMS_meth_pheno = getG2methAtCpGs(EffectsDF_ANNOT$DMS)

## PCA of the methylation
resPCA.allDMS = getResPCA(DMS_allDMS_meth_pheno)

## Run model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
resMod_PCA.allDMS = getModPCA(resPCA.allDMS)

##############
# Model found:
#   BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex) + PCA1:PAT

#     Eliminated  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)   
# PCA1:PAT    0 23868.7 23868.7     1   107  7.3048 0.0080 **

## Scree plot
plotScreePlot(resPCA.allDMS)
## Mainly 3 first component, 4,5 a bit, then big drop

## PCA plot
fig4a <- myPlotPCA(resPCA.allDMS, resMod_PCA.allDMS) + 
  ggtitle("Individuals PCA based on all DMS")
fig4a 

### How much of the BCI variance is explained by each variables?
#   BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex) + PCA1:PAT
modFULL = lmer(BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 PCA1:PAT, data = resMod_PCA.allDMS$metadata)
message("R2c= conditional R2 value associated with fixed effects plus the random effects.")

mod_noPAT = lmer(BCI ~ PCA1 + (1 | brotherPairID) + (1 | Sex),
                 data = resMod_PCA.allDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPAT)[2])*100,2),
  "% of the variance in associated with the paternal infection"))

mod_noPCA1 = lmer(BCI ~ PAT + (1 | brotherPairID) + (1 | Sex),
                  data = resMod_PCA.allDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPCA1)[2])*100,2),
  "% of the variance in associated with the first PCA axis"))

### Plot of the model 
plotModBCI_PCA_allDMS = plotModBCI_PCA(modFULL, resMod_PCA.allDMS, Signifterms = c("PCA1", "PAT"))

plotModBCI_PCA_allDMS$p1
fig4d <- plotModBCI_PCA_allDMS$p2+ 
  ggtitle("Predicted and observed BCI based on PCA1 at all DMS")
fig4d

#########################
## ONLY intergenerational
### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values
DMS_intergenerationalDMS_meth_pheno = getG2methAtCpGs(EffectsDF_ANNOT$DMS[
  EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL"])

## PCA of the methylation
resPCA.intergenerationalDMS = getResPCA(DMS_intergenerationalDMS_meth_pheno)

## Run model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
resMod_PCA.intergenerationalDMS = getModPCA(resPCA.intergenerationalDMS)

##############
# Model found:
#   BCI ~ PCA2 + PAT + (1 | brotherPairID) + (1 | Sex) + PCA2:PAT

#     Eliminated  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)   
# PCA2:PAT     0 29923.5 29923.5     1   107  9.5262 0.002579 **
  
## Scree plot
plotScreePlot(resPCA.intergenerationalDMS)
## Mainly 3 first component, 4,5 a bit, then big drop

## PCA plot
fig4b <- myPlotPCA(resPCA.intergenerationalDMS, resMod_PCA.intergenerationalDMS)+ 
  ggtitle("Individuals PCA based on intergenerational DMS")
fig4b

### How much of the BCI variance is explained by each variables?
modFULL = lmer(BCI ~ PCA2 + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 PCA2:PAT, data = resMod_PCA.intergenerationalDMS$metadata)
message("R2c= conditional R2 value associated with fixed effects plus the random effects.")

mod_noPAT = lmer(BCI ~ PCA2 + (1 | brotherPairID) + (1 | Sex),
                 data = resMod_PCA.intergenerationalDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPAT)[2])*100,2),
  "% of the variance in associated with the paternal infection"))

mod_noPCA2 = lmer(BCI ~ PAT + (1 | brotherPairID) + (1 | Sex),
                  data = resMod_PCA.intergenerationalDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPCA1)[2])*100,2),
  "% of the variance in associated with the second PCA axis"))

### Plot of the model 
plotModBCI_PCA_intergenerationalDMS = plotModBCI_PCA(
  modFULL, resMod_PCA.intergenerationalDMS, Signifterms = c("PCA2", "PAT"))

plotModBCI_PCA_intergenerationalDMS$p1
fig4e <- plotModBCI_PCA_intergenerationalDMS$p2 + 
  ggtitle("Predicted and observed BCI based on PCA2 at intergenerational DMS")
fig4e

#########################
## ONLY infection-induced
### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values
DMS_infectionInducedDMS_meth_pheno = getG2methAtCpGs(EffectsDF_ANNOT$DMS[
  EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED"])

## PCA of the methylation
resPCA.infectionInducedDMS = getResPCA(DMS_infectionInducedDMS_meth_pheno)

## Run model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
resMod_PCA.infectionInducedDMS = getModPCA(resPCA.infectionInducedDMS)

##############
# Model found:
#   BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex) 

#     Eliminated  Sum Sq Mean Sq NumDF DenDF F value Pr(>F)   
# PCA1                            0  58946   58946     1   108 17.4717 5.943e-05 ***
# PAT                             0  42725   42725     1   108 12.6635 0.0005554 ***

## Scree plot
plotScreePlot(resPCA.infectionInducedDMS)
## Mainly 3 first component, 4,5 a bit, then big drop

## PCA plot
fig4c <- myPlotPCA(resPCA.infectionInducedDMS, resMod_PCA.infectionInducedDMS)+ 
  ggtitle("Individuals PCA based on infection-induced DMS")
fig4c

### How much of the BCI variance is explained by each variables?
modFULL = lmer(BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex), 
                 data = resMod_PCA.infectionInducedDMS$metadata)
message("R2c= conditional R2 value associated with fixed effects plus the random effects.")

mod_noPAT = lmer(BCI ~ PCA1 + (1 | brotherPairID) + (1 | Sex),
                 data = resMod_PCA.infectionInducedDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPAT)[2])*100,2),
  "% of the variance in associated with the paternal infection"))

mod_noPCA1 = lmer(BCI ~ PAT + (1 | brotherPairID) + (1 | Sex),
                  data = resMod_PCA.infectionInducedDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPCA1)[2])*100,2),
  "% of the variance in associated with the first PCA axis"))

### Plot of the model 
plotModBCI_PCA_infectionInducedDMS = plotModBCI_PCA(
  modFULL, resMod_PCA.infectionInducedDMS, Signifterms = c("PCA1", "PAT"))

plotModBCI_PCA_infectionInducedDMS$p1
fig4f <- plotModBCI_PCA_infectionInducedDMS$p2  + 
  ggtitle("Predicted and observed BCI based on PCA1 at infection-induced DMS")
fig4f

# save
pdf(file = "../../dataOut/fig/Fig4_phenoMethPlot.pdf", width = 15, height = 10)
gridExtra::grid.arrange(fig4a + theme(legend.position = c(0.8, 0.2)),
                        fig4b + theme(legend.position = "none"),
                        fig4c + theme(legend.position = "none"),
                        fig4d + guides(fill = "none", linetype = guide_legend()),
                        fig4e + guides(fill = "none", linetype = guide_legend()),
                        fig4f + guides(fill = "none", linetype = guide_legend()),
                        ncol=3)
dev.off()

###################################
## Check the most interesting genes

# The function dimdesc() can be used to identify the most correlated variables with a given principal component.
mydimdesc = dimdesc(myPCA_DMS_allDMS$res.PCA, axes = c(1,2), proba = 0.05)

print(paste(nrow(mydimdesc$Dim.1$quanti), "CpG sites most correlated (p < 0.05) with the first principal component"))
print(paste(nrow(mydimdesc$Dim.2$quanti), "CpG sites most correlated (p < 0.05) with the second principal component"))

## The most correlated with PCA1:
EffectsDF_ANNOT[EffectsDF_ANNOT$DMS %in% rownames(mydimdesc$Dim.1$quanti),] %>% 
  dplyr::select(c("GeneSymbol", "feature.name", "Note", "chrom", "nDMSperGenekb", 
                  "GeneName", "Function", "effect"))%>%
  unique -> annotPCA1

annotPCA1$feature.name %>% length #482 genes

write.csv(annotPCA1, "../../dataOut/fig/TableS2_annotPCA1_482genes.csv", row.names = F)

### GO term for these CpGs
dfGO_PCA1_482 = makedfGO(
  annot = EffectsDF_ANNOT[EffectsDF_ANNOT$DMS %in% rownames(mydimdesc$Dim.1$quanti),] %>%
    distinct(feature.name,.keep_all = TRUE), 
  gene_universe = gene_universe,
  effect = "PCAaxis1",
  label = "PCAaxis1")

pdf(file = "../../dataOut/fig/FigS6_GOplot_482genesPCA2.pdf", width = 30, height = 4)
makeGOplot(dfGO_PCA1_482)
dev.off()

## With slim GO terms
dfGO_PCA1_482_SLIM = makeGOslim(dfGO = dfGO_PCA1_482)

## Find genes that seemed important
EffectsDF_ANNOT[EffectsDF_ANNOT$DMS %in% rownames(mydimdesc$Dim.1$quanti) &
                  EffectsDF_ANNOT$GeneSymbol %in% "PPFIBP1",]

message("R11 done. \n")

