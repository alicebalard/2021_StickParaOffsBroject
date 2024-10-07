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

## Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
## PCA of the methylation
myPCA_DMS_allDMS = DMS_allDMS_meth_pheno %>% `rownames<-`(.[,1]) %>% 
  dplyr::select(matches("Gy_chr")) %>% 
  as.matrix %>% myPCA_mod(incomplete = T)

# Model found:
# BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + PCA1:No.Worms + No.Worms:PAT
#                         Eliminated  Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
# PCA1:No.Worms                   0 13362.9 13362.9     1   105  4.0190 0.047563 * 
# No.Worms:PAT                    0 23271.6 23271.6     1   105  6.9992 0.009408 **

### How much of the BCI variance is explained by each variables?
# BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + PCA1:No.Worms + No.Worms:PAT
modFULL = lmer(BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 No.Worms:PCA1 + No.Worms:PAT, data = myPCA_DMS_allDMS$metadata)

message("R2c= conditional R2 value associated with fixed effects plus the random effects.")

mod_noPAT = lmer(BCI ~ PCA1 + No.Worms + (1 | brotherPairID) + (1 | Sex) + 
                   No.Worms:PCA1, data = myPCA_DMS_allDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPAT)[2])*100,2),
  "% of the variance in associated with the paternal infection"))

mod_noWorms = lmer(BCI ~ PCA1 + PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = myPCA_DMS_allDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noWorms)[2])*100,2),
  "% of the variance in associated with the number of worms"))

mod_noPCA1 = lmer(BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                    No.Worms:PAT, data = myPCA_DMS_allDMS$metadata)
message(paste0(
  round((MuMIn::r.squaredGLMM(modFULL)[2] - MuMIn::r.squaredGLMM(mod_noPCA1)[2])*100,2),
  "% of the variance in associated with the first PCA axis"))

### Plot of the model
# No.Worms:PAT
p1 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PAT"), alpha=.4)+ 
  scale_color_manual(values = setNames(colOffs, NULL)[1:2])+
  scale_fill_manual(values = setNames(colOffs, NULL)[1:2])+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:paternal treatment interaction")

# No.Worms:PAT
p2 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PCA1 [-10, 10]"))+ 
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("black","red"))+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:PCA1 interaction")

# save
pdf(file = "../../dataOut/fig/Fig4_phenoMethPlot_alleffects.pdf", width = 8, height = 5)
gridExtra::grid.arrange(p1,p2, ncol=2)
dev.off()

#########################
## ONLY intergenerational
## PCA1:worms close to significant!

### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values
DMS_intergenerational_meth_pheno = getG2methAtCpGs(EffectsDF_ANNOT$DMS[
  EffectsDF_ANNOT$effect %in% "INTERGENERATIONAL"])

## Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
## PCA of the methylation
myPCA_DMS_intergenerational = DMS_intergenerational_meth_pheno %>% `rownames<-`(.[,1]) %>% 
  dplyr::select(matches("Gy_chr")) %>% 
  as.matrix %>% myPCA_mod(incomplete = T)

# Model found:
# BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
# PCA1:No.Worms                  11 10636.8 10636.8     1   105  3.1309 0.07972 .

## NB: not the model selected! but to compare with the overall
modFULL = lmer(BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 No.Worms:PCA1 + No.Worms:PAT, data = myPCA_DMS_intergenerational$metadata)

### Plot of the model
# No.Worms:PAT
p1 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PAT"), alpha=.4)+ 
  scale_color_manual(values = setNames(colOffs, NULL)[1:2])+
  scale_fill_manual(values = setNames(colOffs, NULL)[1:2])+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:paternal treatment interaction")

# No.Worms:PAT
p2 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PCA1 [-10, 10]"))+ 
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("black","red"))+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:PCA1 interaction")

# save
pdf(file = "../../dataOut/fig/FigS4B_phenoMethPlot_intergenerational.pdf", width = 8, height = 5)
gridExtra::grid.arrange(p1,p2, ncol=2)
dev.off()

#########################
## ONLY infection-induced 
## PCA1:worms very far from significant!

### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values
DMS_infectioninduced_meth_pheno = getG2methAtCpGs(EffectsDF_ANNOT$DMS[
  EffectsDF_ANNOT$effect %in% "INFECTION_INDUCED"])

## Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
## PCA of the methylation
myPCA_DMS_infectioninduced = DMS_infectioninduced_meth_pheno %>% `rownames<-`(.[,1]) %>% 
  dplyr::select(matches("Gy_chr")) %>% 
  as.matrix %>% myPCA_mod(incomplete = T)

# Model found:
# BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
# PCA1:No.Worms                   8  3346.5  3346.5     1   102  0.9859 0.32311  

## NB: not the model selected! but to compare with the overall
modFULL = lmer(BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 No.Worms:PCA1 + No.Worms:PAT, data = myPCA_DMS_infectioninduced$metadata)

### Plot of the model
# No.Worms:PAT
p1 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PAT"), alpha=.4)+ 
  scale_color_manual(values = setNames(colOffs, NULL)[1:2])+
  scale_fill_manual(values = setNames(colOffs, NULL)[1:2])+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:paternal treatment interaction")

# No.Worms:PAT
p2 = plot_model(modFULL, type = "pred", terms = c("No.Worms","PCA1 [-10, 10]"))+ 
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("black","red"))+
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of BCI", subtitle = "No.Worms:PCA1 interaction")

# save
pdf(file = "../../dataOut/fig/FigS4A_phenoMethPlot_infectioninduced.pdf", width = 8, height = 5)
gridExtra::grid.arrange(p1,p2, ncol=2)
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
