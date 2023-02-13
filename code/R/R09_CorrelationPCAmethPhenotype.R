# Each script sources the previous script of the pipeline if needed
source("R08_DMSannotation.R")

## Correlation between methylation (after PCA) and phenotype (nbr worms, BCI)

# Approach:

#  For each DMS G1 and/or G2 effects:
# 1. Extract methylation values: raw beta values at DMS shared by \>4 (or more) BP
# 2. PCA
# 3. Extract axes 1 & 2
# 4. Correlation with parasite load/BCI
# 5. If result positive, annotate the CpG associated with significant axis

### PCA based on all methylation values at DMS positions detected for all effects, with imputation of missing values

# Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
RESPCA <- getPCACpG(DMSvec=unique(c(DMS_PaternalEffect_4BPmin, DMS_OffspringEffect_4BPmin)), effect="all effects")
## Save for later 
RESPCAgeneral <- RESPCA
# 2273 DMS linked with all effects
# [1] "The chosen model is:"
# BCI ~ PCA2 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
#     PCA2:PAT + No.Worms:PAT
# PCA2:PAT                        0  40035   40035     1   105 12.6251 0.0005717 ***
# No.Worms:PAT                    0  19851   19851     1   105  6.2601 0.0138920 *  
# [1] "1169 CpG sites most correlated (p < 0.05) with the first principal component"
# [1] "1137 CpG sites most correlated (p < 0.05) with the second principal component"

formula(RESPCA$PCA_percAtDMS_imputed$modSel)
# The SECOND PCA axis is significant in BCI
# [1] "The chosen model is:"
# BCI ~ PCA2 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
#   PCA2:PAT + No.Worms:PAT

### How much of the BCI variance is explained by each variables?
mod_noworms = lmer(BCI ~ PCA2 + PAT + PCA2:PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPAT = lmer(BCI ~ PCA2 + No.Worms + (1 | brotherPairID) + (1 | Sex), 
                 data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPCA2 = lmer(BCI ~ No.Worms + PAT + No.Worms:PAT +(1 | brotherPairID) + (1 | Sex), 
                  data = RESPCA$PCA_percAtDMS_imputed$metadata)

# R2c conditional R2 value associated with fixed effects plus the random effects.
A = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noworms)[2])*100
B = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
C = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPCA2)[2])*100
round(A, 2) #10.72% of the variance in associated with the parasite load (number of worms)
round(B, 2) #21.99% of the variance in associated with the paternal treatment
round(C, 2) #9.46% of the variance in associated with the second PCA axis

### Plot of the model
phenoMethPlot <- plot(ggpredict(RESPCA$PCA_percAtDMS_imputed$modSel, terms = c("No.Worms", "PCA2", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring")
phenoMethPlot

# save
pdf(file = "../../dataOut/phenotypeMeth/phenoMethPlot_alleffects.pdf", width = 7, height = 5)
phenoMethPlot
dev.off()

########################
### PCA based on all methylation values at DMS positions detected for each effects separately, with imputation of missing values

# Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
RESPCA <- getPCACpG(DMSvec=DMS_G1onlyEffect_4BPmin, effect="G1")
formula(RESPCA$PCA_percAtDMS_imputed$modSel)

### How much of the BCI variance is explained by each variables?
mod_noworms = lmer(BCI ~ PCA2 + PAT + PCA2:PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPAT = lmer(BCI ~ PCA2 + No.Worms + (1 | brotherPairID) + (1 | Sex), 
                 data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPCA2 = lmer(BCI ~ + No.Worms + PAT + No.Worms:PAT +(1 | brotherPairID) + (1 | Sex), 
                  data = RESPCA$PCA_percAtDMS_imputed$metadata)

# R2c conditional R2 value associated with fixed effects plus the random effects.
A = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noworms)[2])*100
B = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
C = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPCA2)[2])*100
round(A, 2) #10.41% of the variance in associated with the parasite load (number of worms)
round(B, 2) #22.02% of the variance in associated with the paternal treatment
round(C, 2) #9.6% of the variance in associated with the second PCA axis

## Summary:
# 1640 DMS linked with intergenerational effect
# 862 CpG sites most correlated (p < 0.05) with PCA1
# 826 CpG sites most correlated (p < 0.05) with PCA2
# The chosen model is:
# BCI ~ PCA2 + No.Worms + paternal treatment + PCA2:paternal treatment + No.Worms:paternal treatment + (1 | father's family) + (1 | Sex)
# Backward reduced elimination with Satterthwaite method
# PCA2:paternal treatment F=12.48, p-value<0.001
# No.Worms:paternal treatment F=6.47, p-value=0.012
# Variance in BCI associated with:
# No.Worms=10.4% paternal treatment=22% PCA2=9.6%

### Plot of the model
phenoMethPlotG1 <- plot(ggpredict(RESPCA$PCA_percAtDMS_imputed$modSel, terms = c("No.Worms", "PCA2", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring: G1")
phenoMethPlotG1

# save
pdf(file = "../../dataOut/phenotypeMeth/phenoMethPlot_G1.pdf", width = 7, height = 5)
phenoMethPlotG1
dev.off()

# Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
RESPCA <- getPCACpG(DMSvec=DMS_G2onlyEffect_4BPmin, effect="G2")

formula(RESPCA$PCA_percAtDMS_imputed$modSel)

### How much of the BCI variance is explained by each variables?
mod_noworms = lmer(BCI ~ PCA2 + PAT + PCA2:PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPAT = lmer(BCI ~ PCA2 + No.Worms + (1 | brotherPairID) + (1 | Sex), 
                 data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPCA2 = lmer(BCI ~ + No.Worms + PAT + No.Worms:PAT +(1 | brotherPairID) + (1 | Sex), 
                  data = RESPCA$PCA_percAtDMS_imputed$metadata)

# R2c conditional R2 value associated with fixed effects plus the random effects.
A = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noworms)[2])*100
B = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
C = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPCA2)[2])*100
round(A, 2) #11.81% of the variance in associated with the parasite load (number of worms)
round(B, 2) #18.88% of the variance in associated with the paternal treatment
round(C, 2) #6.08% of the variance in associated with the second PCA axis

## Summary:
# 309 DMS linked with infection-induced effect
# 181 CpG sites most correlated (p < 0.05) with PCA1
# 178 CpG sites most correlated (p < 0.05) with PCA2
# The chosen model is:
# BCI ~ PCA2 + No.Worms + paternal treatment + PCA2:paternal treatment + No.Worms:paternal treatment + (1 | father's family) + (1 | Sex)
# Backward reduced elimination with Satterthwaite method
# PCA2:paternal treatment F=12.48, p-value<0.001
# No.Worms:paternal treatment F=6.47, p-value=0.012
# Variance in BCI associated with:
# No.Worms=11.8% paternal treatment=18.9% PCA2=6.1%

### Plot of the model
phenoMethPlotG2 <- plot(ggpredict(RESPCA$PCA_percAtDMS_imputed$modSel, terms = c("No.Worms", "PCA2", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring: G2")
phenoMethPlotG2

# save
pdf(file = "../../dataOut/phenotypeMeth/phenoMethPlot_G2.pdf", width = 7, height = 5)
phenoMethPlotG2
dev.off()

# Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
RESPCA <- getPCACpG(DMSvec=DMS_G1G2additiveEffect_4BPmin, effect="additive")

formula(RESPCA$PCA_percAtDMS_imputed$modSel)

### How much of the BCI variance is explained by each variables?
mod_noworms = lmer(BCI ~ PCA2 + PAT + PCA2:PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPAT = lmer(BCI ~ PCA2 + No.Worms + (1 | brotherPairID) + (1 | Sex), 
                 data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPCA2 = lmer(BCI ~ + No.Worms + PAT + No.Worms:PAT +(1 | brotherPairID) + (1 | Sex), 
                  data = RESPCA$PCA_percAtDMS_imputed$metadata)

# R2c conditional R2 value associated with fixed effects plus the random effects.
A = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noworms)[2])*100
B = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
C = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPCA2)[2])*100
round(A, 2) #10.86% of the variance in associated with the parasite load (number of worms)
round(B, 2) #20.03% of the variance in associated with the paternal treatment
round(C, 2) #9.17% of the variance in associated with the second PCA axis

## Summary:
# 173 DMS linked with additive effect
# 95 CpG sites most correlated (p < 0.05) with PCA1
# 81 CpG sites most correlated (p < 0.05) with PCA2
# The chosen model is:
# BCI ~ PCA2 + No.Worms + paternal treatment + PCA2:paternal treatment + No.Worms:paternal treatment + (1 | father's family) + (1 | Sex)
# Backward reduced elimination with Satterthwaite method
# PCA2:paternal treatment F=12.48, p-value<0.001
# No.Worms:paternal treatment F=6.47, p-value=0.012
# Variance in BCI associated with:
# No.Worms=10.9% paternal treatment=20% PCA2=9.2%

### Plot of the model
phenoMethPlotadditive <- plot(ggpredict(RESPCA$PCA_percAtDMS_imputed$modSel, terms = c("No.Worms", "PCA2", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring: additive")
phenoMethPlotadditive

# save
pdf(file = "../../dataOut/phenotypeMeth/phenoMethPlot_additive.pdf", width = 7, height = 5)
phenoMethPlotadditive
dev.off()

# Make PCA and model lmer(BCI ~ PCA1*PCA2*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=metadata)
RESPCA <- getPCACpG(DMSvec=DMS_G1G2interactionEffect_4BPmin, effect="interaction")

formula(RESPCA$PCA_percAtDMS_imputed$modSel)

### How much of the BCI variance is explained by each variables?
mod_noworms = lmer(BCI ~ PCA2 + PAT + PCA2:PAT + (1 | brotherPairID) + (1 | Sex), 
                   data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPAT = lmer(BCI ~ PCA2 + No.Worms +  PCA2:No.Worms +(1 | brotherPairID) + (1 | Sex), 
                 data = RESPCA$PCA_percAtDMS_imputed$metadata)
mod_noPCA2 = lmer(BCI ~ No.Worms + PAT + No.Worms:PAT +(1 | brotherPairID) + (1 | Sex), 
                  data = RESPCA$PCA_percAtDMS_imputed$metadata)

# R2c conditional R2 value associated with fixed effects plus the random effects.
A = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noworms)[2])*100
B = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
C = (MuMIn::r.squaredGLMM(RESPCA$PCA_percAtDMS_imputed$modSel)[2] -
       MuMIn::r.squaredGLMM(mod_noPCA2)[2])*100
round(A, 2) #14.75% of the variance in associated with the parasite load (number of worms)
round(B, 2) #19.02% of the variance in associated with the paternal treatment
round(C, 2) #7.05% of the variance in associated with the second PCA axis

## Summary:
# 151 DMS linked with interaction effect
# 72 CpG sites most correlated (p < 0.05) with PCA1
# 54 CpG sites most correlated (p < 0.05) with PCA2
# The chosen model is:
# BCI ~ PCA2 + No.Worms + paternal treatment +
# PCA2:No.Worms + PCA2:paternal treatment + No.Worms:paternal treatment + PCA2:No.Worms:paternal treatment +
# (1 | father's family) + (1 | Sex)
# Backward reduced elimination with Satterthwaite method
# PCA2:paternal treatment F=12.48, p-value<0.001
# No.Worms:paternal treatment F=6.47, p-value=0.012
# Variance in BCI associated with:
# No.Worms=14.8% paternal treatment=19% PCA2=7.1%

### Plot of the model
phenoMethPlotinteraction <- plot(ggpredict(RESPCA$PCA_percAtDMS_imputed$modSel, terms = c("No.Worms", "PCA2", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring: interaction")
phenoMethPlotinteraction

# save
pdf(file = "../../dataOut/phenotypeMeth/phenoMethPlot_interaction.pdf", width = 7, height = 5)
phenoMethPlotinteraction
dev.off()

################
### Annotate the genes linked with axis 2 of the PCA
annotPCAaxisFull <- myHomebrewDMSannotation(DMSvec = paste(RESPCAgeneral$CpGPCA2$chr, RESPCAgeneral$CpGPCA2$end),
                                            myannotBed12 = annotBed12, myannotGff3 = annotGff3)

# merge with full table to add effect
annotPCAaxisFull = merge(annotPCAaxisFull, allDMSAnnot)

annotPCAaxis = annotPCAaxisFull %>% 
  dplyr::select(c("GeneSymbol", "feature.name", "Note", "chrom", "nDMSperGenekb", "ENTREZID", "description", "summary", "effect"))%>% 
  unique

write.csv(annotPCAaxis, "../../dataOut/annotPCA2_1137DMS_437genes_supTabS2.csv", row.names = F)

### Plot annotated Manhattan plots for those 437 genes associated with PCA2!
P=plotManhattanGenesDMS(annotFile = annotPCAaxisFull, GYgynogff = GYgynogff)
P

pdf("../../dataOut/ManhattanPlots1137DMS437genes.pdf", width = 10, height = 3)
P
dev.off()

### GO term for these CpGs
GO_PCA2_1137DMS = makedfGO(annotPCAaxisFull, gene_universe, effect = "PCAaxis2")

GOplot <- GO_PCA2_1137DMS %>% ggplot(aes(x=Effect, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name="adjusted\np-value", low = "red", high = "blue") +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), legend.position="left") + # grey box for legend
  facet_grid(fct_inorder(GO.category)~., scales="free",space = "free")+
  scale_y_discrete(limits=rev) # revers axis to have alphabetical order

GOplot

pdf(GOplot, file = "../../dataOut/GOplotPCA2_supF6.pdf", width = 6, height = 8)
GOplot
dev.off()

################
# GO of 437 genes linked with the 1137 DMS PCA2, split by effect (NB some genes in different effects; genes with at least 1 CpG in one effect)
GO_G1_sub = makedfGO(annotPCAaxisFull[!is.na(annotPCAaxisFull$G1),], gene_universe, effect = "G1")
GO_G2_sub = makedfGO(annotPCAaxisFull[!is.na(annotPCAaxisFull$G2),], gene_universe, effect = "G2")
GO_addit_sub = makedfGO(annotPCAaxisFull[!is.na(annotPCAaxisFull$addit),], gene_universe, effect = "addit")
GO_inter_sub = makedfGO(annotPCAaxisFull[!is.na(annotPCAaxisFull$inter),], gene_universe, effect = "inter")

dfGO_sub = rbind(GO_G1_sub, GO_G2_sub, GO_addit_sub, GO_inter_sub)

### GO plot
GOplot_sub <- dfGO_sub %>% ggplot(aes(x=Effect, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name="adjusted\np-value", low = "red", high = "blue") +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("Treatments comparison") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), legend.position="left") + # grey box for legend
  facet_grid(fct_inorder(GO.category)~., scales="free",space = "free")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) # split too long GO names in half
GOplot_sub

pdf(GOplot_sub, file = "../../dataOut/GOplot4Venncat_split.pdf", width = 6, height = 15)
GOplot_sub
dev.off()

#############################
## Test of specific positions
### First, some plots

# Top 5 genes with more DMS
# G1: DMS on MTZ2, HMX1, GTF2IRD2, PHACTR3, SPO11 -> on PCA2 only HMX1 (2), GTF2IRD2 (1), SPO11 (2)
# G2: DMS on TRIM16, ZNF691, DCT, ACAA1, MGAT2 -> on PCA2 only TRIM16 (8), DCT (5), ACAA1 (2), MGAT2 (1)
# Additive: ZNF518B, ERI2, GYAR, ZMYND19, SLC7A6 -> on PCA2 only GYAR (1), SLC7A6 (1)
# Interaction: CIAO2B, CYB561D2, BCL7A, SNIP1, TENT5A -> on PCA2 only CIAO2B (1), CYB561D2 (1), BCL7A (1), 
# To find out: annotPCAaxisFull[annotPCAaxisFull$GeneSymbol %in% "TENT5A",]

mygene = "GYAR"

mypos = paste(annotPCAaxisFull[annotPCAaxisFull$GeneSymbol %in% mygene,"chrom"],
              annotPCAaxisFull[annotPCAaxisFull$GeneSymbol %in% mygene,"start"])

mydf = methylKit::select(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
                         which(paste(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr, uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$start) %in% 
                                 mypos[1])) %>%
  methylKit::percMethylation() %>% melt %>% dplyr::select(c("Var2", "value")) %>% dplyr::rename("SampleID"="Var2") %>%
  merge(fullMetadata_OFFS)

mydf = na.omit(mydf[c("value", "BCI", "No.Worms", "PAT", "brotherPairID", "Sex")]) 

mymod = lmerTest::lmer(BCI ~ value*No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=mydf)

mymodSel=lmer(formula = attr(attr(lmerTest::step(mymod, reduce.random = F), "drop1"), "heading")[3],
     data=mydf, REML = F)

mymodSel%>%formula

plot(ggpredict(mymodSel, terms = c("No.Worms","value", "PAT")), add.data = TRUE, alpha = .08) +
  theme_bw() +
  scale_color_gradient(low = "white", high = "red")+
  scale_fill_gradient(low = "white", high = "red") +
  ylab("Body Condition Index") + xlab("Number of worms")+
  ggtitle("Predicted values of Body Condition Index in offspring")
