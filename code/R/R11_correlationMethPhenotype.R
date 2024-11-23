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

############
## Print PCA
barplot(myPCA_DMS_allDMS$res.PCA$eig[, 2], names.arg=1:nrow(myPCA_DMS_allDMS$res.PCA$eig), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(myPCA_DMS_allDMS$res.PCA$eig), myPCA_DMS_allDMS$res.PCA$eig[, 2], 
      type="b", pch=19, col = "red")
## Mainly 3 first component, 4,5 a bit, then big drop

pcaplot0 <-fviz_pca_ind(myPCA_DMS_allDMS$res.PCA, label="none", 
             habillage=fullMetadata_OFFS$trtG1G2[
               match(rownames(myPCA_DMS_allDMS$res.PCA$ind$coord), fullMetadata_OFFS$SampleID)], 
             pointsize =3, addEllipses=TRUE)+
  scale_color_manual(values = colOffs)+
  scale_fill_manual(values = colOffs)
pcaplot0

##############
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
makeplotModel <- function(modFULL, mymin){
  # No.Worms:PAT
  pred_data <- ggpredict(modFULL, terms = c("No.Worms", "PAT"))
  
  dfobs = myPCA_DMS_allDMS$metadata
  dfobs$x = dfobs$No.Worms
  dfobs$predicted = dfobs$BCI
  dfobs$group = dfobs$PAT
  
  p1 = ggplot(pred_data, aes(x = x, y = predicted, fill = group, color = group)) +
    geom_line(aes(group = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
                alpha = 0.2, linetype = "blank")+
    geom_jitter(data=dfobs, aes(fill = group), size = 3,
                height = 0, width = .1)+
    labs(x = "Number of Worms", y = "Body Condition Index", 
         color = "PAT", fill = "PAT") +
    ggtitle("Predicted and observed values of BCI", 
            subtitle = "No.Worms:paternal treatment interaction")+
    scale_color_manual(values = setNames(colOffs, NULL)[1:2])+
    scale_fill_manual(values = setNames(colOffs, NULL)[1:2])+
    ylim(mymin,300)  # Set y-axis limits
  
  # No.Worms:PCA1
  pred_data <- ggpredict(modFULL, terms = c("No.Worms", "PCA1 [-18, 13]"))
  pred_data$group <- as.numeric(as.character(pred_data$group))
  
  dfobs = myPCA_DMS_allDMS$metadata
  dfobs$x = dfobs$No.Worms
  dfobs$predicted = dfobs$BCI
  dfobs$group = dfobs$PCA1
  
  p2 = ggplot(pred_data, aes(x = x, y = predicted, fill = group, color = group)) +
    geom_line(aes(group = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
                alpha = 0.2, linetype = "blank")+
    geom_jitter(data=dfobs, aes(fill = group), size = 3,
                height = 0, width = .1)+
    labs(x = "Number of Worms", y = "Body Condition Index", 
         color = "PCA1", fill = "PCA1") +
    ggtitle("Predicted and observed values of BCI", 
            subtitle = "No.Worms:PCA1 interaction")+
    scale_color_gradient(low = "black", high = "red")+
    scale_fill_gradient(low = "black", high = "red")+
    ylim(mymin,300)  # Set y-axis limits
  
  return(list(p1=p1, p2=p2))
}

p1 = makeplotModel(modFULL, mymin=-380)$p1
p2 = makeplotModel(modFULL, mymin=-380)$p2

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

############
## Print PCA
barplot(myPCA_DMS_intergenerational$res.PCA$eig[, 2], 
        names.arg=1:nrow(myPCA_DMS_intergenerational$res.PCA$eig), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(myPCA_DMS_intergenerational$res.PCA$eig), 
      myPCA_DMS_intergenerational$res.PCA$eig[, 2], 
      type="b", pch=19, col = "red")
## Mainly 2 first components, 3,4 a bit, then big drop

pcaplot1 <- fviz_pca_ind(myPCA_DMS_intergenerational$res.PCA, label="none", 
             habillage=fullMetadata_OFFS$trtG1G2[
               match(rownames(myPCA_DMS_intergenerational$res.PCA$ind$coord), 
                     fullMetadata_OFFS$SampleID)], 
             pointsize =3, addEllipses=TRUE)+
  scale_color_manual(values = colOffs)+
  scale_fill_manual(values = colOffs)+
  ggtitle("Individuals - PCA of intergenerational DMS")
pcaplot1

#############
# Model found:
# BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
# PCA1:No.Worms                  11 10636.8 10636.8     1   105  3.1309 0.07972 .

## NB: not the model selected! but to compare with the overall
modFULL = lmer(BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 No.Worms:PCA1 + No.Worms:PAT, data = myPCA_DMS_intergenerational$metadata)

### Plot of the model
p1 = makeplotModel(modFULL, mymin=-500)$p1
p2 = makeplotModel(modFULL,mymin=-500)$p2

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

############
## Print PCA
barplot(myPCA_DMS_infectioninduced$res.PCA$eig[, 2], 
        names.arg=1:nrow(myPCA_DMS_infectioninduced$res.PCA$eig), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(myPCA_DMS_infectioninduced$res.PCA$eig), 
      myPCA_DMS_infectioninduced$res.PCA$eig[, 2], 
      type="b", pch=19, col = "red")
## Mainly 2 first components, 3,4,5 a bit, then big drop

pcaplot2 <- fviz_pca_ind(myPCA_DMS_infectioninduced$res.PCA, label="none", 
             habillage=fullMetadata_OFFS$trtG1G2[
               match(rownames(myPCA_DMS_infectioninduced$res.PCA$ind$coord), 
                     fullMetadata_OFFS$SampleID)], 
             pointsize =3, addEllipses=TRUE)+
  scale_color_manual(values = colOffs)+
  scale_fill_manual(values = colOffs)+
  ggtitle("Individuals - PCA of infection-induced DMS")
pcaplot2

#############
# Model found:
# BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
# PCA1:No.Worms                   8  3346.5  3346.5     1   102  0.9859 0.32311  

## NB: not the model selected! but to compare with the overall
modFULL = lmer(BCI ~ PCA1 + No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + 
                 No.Worms:PCA1 + No.Worms:PAT, data = myPCA_DMS_infectioninduced$metadata)

### Plot of the model
p1 = makeplotModel(modFULL, mymin=-500)$p1
p2 = makeplotModel(modFULL,mymin=-500)$p2

# save
pdf(file = "../../dataOut/fig/FigS4A_phenoMethPlot_infectioninduced.pdf", width = 10, height = 5)
gridExtra::grid.arrange(p1,p2, ncol=2)
dev.off()

# save PCA plots
pdf(file = "../../dataOut/fig/Figxx_PCAplots.pdf", width = 10, height = 6)
gridExtra::grid.arrange(pcaplot0,pcaplot1,pcaplot2, ncol=2)
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

