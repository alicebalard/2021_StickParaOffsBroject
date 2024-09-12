## Global methylation and fitness
## A. Balard
## September 2024
## Produces fig S4

# Each script sources the previous script of the pipeline if needed
source("R04_getDiffMeth_PQLseq_runInCLUSTER.R") 

print("Number of CpG positions shared by all fish:")
print(nrow(uniteCovALL_woSexAndUnknowChr))# 60861 

print("Number of CpG positions shared by half the parents in each trt group, overlapping with the offspring ones:")
print(nrow(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP))# 1002565 

print("Number of CpG positions shared by half the offspring in each trt group, overlapping with the paternal ones")
print(nrow(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP)) # 1002565

## Calculate BCI
# Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor of energy reserves and reproductive success, was calculated using there residuals from the regression of body mass on body length (Chellappaet al.1995).
fullMetadata_OFFS$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|brotherPairID), data=fullMetadata_OFFS))

## and for parents (no sex difference, only males):
fullMetadata_PAR$BCI <- residuals(lmer(Wnettofin ~ Slfin + (1|brotherPairID), data=fullMetadata_PAR))

# Effect of paternal treatment on body condition of offspring: Kaufmann et al. 2014: "To investigate in which way paternal G1 exposure affected
# offspring tolerance, we tested how the relationship between G2 body condition and infection intensity was affected by paternal G1 exposure. This was tested in a linear mixed model on G2 body condition with paternal G1 treatment and the interaction between paternal G1 treatment and G2 infection intensity as fixed effects. Maternal half-sibship identity was set as a random effect"

## Effect of paternal exposure on tolerance

### BCI per trt
# Effect of treatment groups of offspring on body condition(Kaufmann et al. 2014): "The linear mixed effect model (nlme function in R) included
# G2 body condition as dependent variable, sex, G2 treatment (exposed vs. control), paternal G1 treatment (exposed vs. control) and their interactions as fixed effects as well as maternal G2 half-sibship identity as a random effect"

mod1 <- lme(BCI ~ offsTrt * patTrt, random=~1|brotherPairID,data=fullMetadata_OFFS)
anova(mod1) # strong significant effect of both offspring trt & paternal + interactions

mod1.2 <- lme(BCI ~  trtG1G2, random=~1|brotherPairID,data=fullMetadata_OFFS)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than
## every other group, same as Kaufmann et al. 2014

myplot1 <- ggplot(fullMetadata_OFFS, aes(x=trtG1G2, y = BCI, fill=trtG1G2))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("NE_control", "NE_exposed")),
              map_signif_level=TRUE, annotations="***",
              y_position = 150, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons = list(c("NE_exposed", "E_control")),
              map_signif_level=TRUE, annotations="***",
              y_position = 200, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons = list(c("NE_exposed", "E_exposed")),
              map_signif_level=TRUE, annotations="***",
              y_position = 250, tip_length = 0, vjust=0.4) +
  scale_fill_manual(values = colOffs)+
  theme_bw() + theme(legend.position = "none") +
  ylab("Body Condition Index") +
  scale_x_discrete(labels=c("NE_control" = "G1 control\nG2 control", "NE_exposed" = "G1 control\nG2 infected",
                            "E_control" = "G1 infected\nG2 control", "E_exposed" = "G1 infected\nG2 infected"),
                   name = NULL)
myplot1

### Tolerance (slope of BCI on worms count) per paternal status
mod_Tol <- lmer(BCI ~ No.Worms*PAT + (1|brotherPairID)+ (1|Sex), data=fullMetadata_OFFS, REML = F)

## Model selection:
step(mod_Tol, reduce.random = F) # Model found: full model

## The slope of BCI on nbrworms varies upon treatment
pdf(file = "../../dataOut/SupplFigS4.pdf", width = 7, height = 5)
plot(ggpredict(mod_Tol, terms = c("No.Worms", "PAT")), show_data=T, jitter = TRUE, 
     dot_alpha =1, alpha = .5)+
  theme_pubr() + # was overridden by plot ggpredict
  ylab("Body Condition Index") + xlab("Number of worms") +
  ggtitle("Predicted values of Body Condition Index in offspring")+
  scale_color_manual(NULL, values = c("#a9d6c1ff", "#d1b000ff")) +
  scale_fill_manual(NULL, values = c("#a9d6c1ff", "#d1b000ff"))  +
  scale_x_continuous(breaks = 0:10)+
  geom_point(size=0)+ # to have color key in legend as point
  guides(colour = guide_legend(override.aes = list(size=3,linetype=0, fill = NA)))
dev.off()

## Link between methylation and fitness (BCI and tolerance)

# Calculate number of methylated sites, mean coverage, and residuals of methylated sites by covered sites (to account for coverage bias)

mycalcRMS <- function(myUniteCov, myMetaData){
  percMethMat = methylKit::percMethylation(myUniteCov)
  # create a dataframe with all info
  percMethDF = data.frame(SampleID = colnames(percMethMat),
                          Nbr_methCpG = colSums(percMethMat>=70 & !is.na(percMethMat)), ## number of methylated sites
                          Nbr_coveredCpG = colSums(!is.na(percMethMat)), ## number of sites covered in this sample
                          Nbr_NOTcoveredCpG = colSums(is.na(percMethMat)),## number of sites NOT covered in this sample
                          MeanCoverage = colMeans(methylKit::getData(myUniteCov)[,myUniteCov@coverage.index], na.rm = T), ## coverage.index: vector denoting which columns in the data correspond to coverage values
                          OverallPercentageMethylation = colMeans(methylKit::percMethylation(myUniteCov), na.rm = T))
  
  ## RMS in this sample based on covered sites
  percMethDF$RMS_coveredCpG = percMethDF$Nbr_methCpG / percMethDF$Nbr_coveredCpG
  ## merge with original metadata:
  myMetaData = merge(myMetaData, percMethDF)
  # calculate also RMS global, considering CpG covered or not (to compare)
  myMetaData$RMS_allCpG_coveredOrNot = myMetaData$Nbr_methCpG / (myMetaData$M.Seqs_rawReads*10e6)
  # calculate residuals of nbr of methCpG by nbr of covered CpG
  myMetaData$res_Nbr_methCpG_Nbr_coveredCpG = residuals(
    lm(myMetaData$Nbr_methCpG ~ myMetaData$Nbr_coveredCpG))
  ## REORDER myMetaData by sample ID
  myMetaData = myMetaData[order(as.numeric(gsub("S", "", myMetaData$SampleID))),]
  return(myMetaData)
}

fullMetadata <- mycalcRMS(uniteCovALL_woSexAndUnknowChr, fullMetadata)

fullMetadata_PAR <- mycalcRMS(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP, fullMetadata_PAR)

fullMetadata_OFFS  <- mycalcRMS(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP, fullMetadata_OFFS)

## Nbr/Ratio of Methylated Sites in different groups
print("mean coverage per CpG site in the full dataset, considering positions covered in all fish")
print(round(mean(fullMetadata$MeanCoverage)))
print("+/-")
print(round(qnorm(0.975)*sd(fullMetadata$MeanCoverage)/sqrt(nrow(fullMetadata)),2))

print("mean coverage per CpG site in G1, considering positions shared by at least half individuals per group, and which overlap with positions retained in G2:")
print(round(mean(fullMetadata_PAR$MeanCoverage)))
print("+/-")
print(round(qnorm(0.975)*sd(fullMetadata_PAR$MeanCoverage)/sqrt(nrow(fullMetadata_PAR)),2))

print("mean coverage per CpG site in G2, considering positions shared by at least half individuals per group, and which overlap with positions retained in G1:")
print(round(mean(fullMetadata_OFFS$MeanCoverage)))
print("+/-")
print(round(qnorm(0.975)*sd(fullMetadata_OFFS$MeanCoverage)/sqrt(nrow(fullMetadata_OFFS)),2))

## Choice of global methylation value
cor.test(fullMetadata_PAR$Nbr_coveredCpG,
         fullMetadata_PAR$Nbr_methCpG, method = "spearman")
## S = 98, p-value <2.2e-06, rho = 0.86

ggplot(fullMetadata_PAR, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

## Check after RMS correction for coverage bias: CORRECTED (p-value = 0.4485)
cor.test(fullMetadata_PAR$Nbr_coveredCpG,
         fullMetadata_PAR$RMS_coveredCpG, method = "spearman")
ggplot(fullMetadata_PAR, aes(x=Nbr_coveredCpG, y=RMS_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

## and with residuals: COMPLETELY CORRECTED p-value = 0.9562
cor.test(fullMetadata_PAR$Nbr_coveredCpG,
         fullMetadata_PAR$res_Nbr_methCpG_Nbr_coveredCpG, method = "spearman")
ggplot(fullMetadata_PAR, aes(x=Nbr_coveredCpG, y=res_Nbr_methCpG_Nbr_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

############
## Offspring:
cor.test(fullMetadata_OFFS$Nbr_coveredCpG,
         fullMetadata_OFFS$Nbr_methCpG, method = "spearman")
ggplot(fullMetadata_OFFS, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  scale_x_continuous("Number of cytosines covered") +
  scale_y_continuous("Number of methylated cytosines") +
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## Plot distance to residuals:
fit <- lm(Nbr_methCpG ~ Nbr_coveredCpG, data = fullMetadata_OFFS)
plotdf <- fullMetadata_OFFS
plotdf$predicted <- predict(fit)   # Save the predicted values
plotdf$residuals <- residuals(fit)
ggplot(plotdf, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_segment(aes(xend = Nbr_coveredCpG, yend = predicted), col = "grey") +
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  scale_x_continuous("Number of cytosines covered") +
  scale_y_continuous("Number of methylated cytosines") +
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## Check after RMS correction for coverage bias: SEMI CORRECTED (p-value = 0.01, rho = -0.24)
cor.test(fullMetadata_OFFS$Nbr_coveredCpG,
         fullMetadata_OFFS$RMS_coveredCpG, method = "spearman")
ggplot(fullMetadata_OFFS, aes(x=Nbr_coveredCpG, y=RMS_coveredCpG))+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  geom_smooth(method = "lm", col="black")+
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## and with residuals: COMPLETELY CORRECTED p-value = 0.51
cor.test(fullMetadata_OFFS$Nbr_coveredCpG,
         fullMetadata_OFFS$res_Nbr_methCpG_Nbr_coveredCpG, method = "spearman")
ggplot(fullMetadata_OFFS, aes(x=Nbr_coveredCpG, y=res_Nbr_methCpG_Nbr_coveredCpG))+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  geom_smooth(method = "lm", col="black")+
  scale_x_continuous("Number of cytosines covered") +
  scale_y_continuous("Residuals of number of methylated cytosines\n on number of cytosines covered") +
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## Why we should we correct for sex

### No difference in mappability (p\>0.05)

mod = lm(`MappingEfficiency%BSBoldvsGynogen` ~ Sex, data = fullMetadata_OFFS)
summary(step(mod))
plot(ggpredict(mod, terms = c("Sex")), add.data = T)
# NB: this is WITH unknown and sex chromosomes, before filtering.

### No difference in number of reads (p\>0.05)

mod = lm(M.Seqs_rawReads ~ Sex, data = fullMetadata_OFFS)
summary(step(mod))
plot(ggpredict(mod, terms = c("Sex")), add.data = T)
# NB: this is WITH unknown and sex chromosomes, before filtering.

### No difference in mean coverage per CpG in the filtered dataset (p\>0.05)
mod = lm(MeanCoverage ~ Sex, data = fullMetadata_OFFS)
summary(step(mod))
plot(ggpredict(mod, terms = c("Sex")), show_data = T)
# NB: this is in G2, considering positions shared by at least 14 fish per treatment group (half individuals per group), and which overlap with positions retained in G1, without sex and unknown chromosome (after filtering) 

### No difference in number of sites covered in the filtered dataset (p\>0.05)
mod = lm(Nbr_coveredCpG ~ Sex, data = fullMetadata_OFFS)
summary(step(mod))
plot(ggpredict(mod, terms = c("Sex")), show_data = T)
# NB: this is in G2, considering positions shared by at least 14 fish per treatment group (half individuals per group), and which overlap with positions retained in G1, without sex and unknown chromosome (after filtering)

### No difference in number of sites covered in the filtered dataset (p\>0.05)
mod = lm(OverallPercentageMethylation ~ Sex, data = fullMetadata_OFFS)
summary(step(mod))
plot(ggpredict(mod, terms = c("Sex")), show_data = T)
#NB: this is in G2, considering positions shared by at least 14 fish per treatment group (half individuals per group), and which overlap with
#positions retained in G1, without sex and unknown chromosome (after filtering)

### Males have a lower global methylation than females (residuals of nbr of methylated sites by nbr of sites covered)
mod = lm(res_Nbr_methCpG_Nbr_coveredCpG ~ Sex, data = fullMetadata_OFFS)
summary(step(mod)) # sex is significant p = 0.000157 ***
anova(mod)

plot(ggpredict(mod, terms = c("Sex")), show_data = T, jitter = T) +
  xlab(NULL)+
  ylab("Residuals of N methylated sites on N covered sites") +
  ggtitle("Predicted values of global methylation in offspring")

## Are mean residuals meth sites different following tolerance slope?

### In all offspring
fullMetadata_OFFS$res_Nbr_methCpG_Nbr_coveredCpG_div1000 <- (fullMetadata_OFFS$res_Nbr_methCpG_Nbr_coveredCpG)/1000

mod_Tol.Meth <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*No.Worms*PAT + (1|brotherPairID)+ (1|Sex),
                     data=fullMetadata_OFFS, REML = F)

## Model selection:
step(mod_Tol.Meth, reduce.random = F) # Model found: BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
## The slope of BCI on nbrworms varies upon treatment but methylation does NOT vary with tolerance
mod_Tol.Meth <- lmer(BCI ~ No.Worms*PAT + (1|brotherPairID)+ (1|Sex),
                     data=fullMetadata_OFFS)

## And by treatment instead of No.worms?
mod_Tol.Meth2 <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*PAT*outcome + (1|brotherPairID)+ (1|Sex),
                      data=fullMetadata_OFFS, REML = F)

## Model selection:
step(mod_Tol.Meth2, reduce.random = F) # Model found: BCI ~ PAT + outcome + (1 | brotherPairID) + (1 | Sex) + PAT:outcome

# The slope of BCI on nbr worms varies upon parental treatment, but methylation does NOT vary with tolerance

### In exposed offspring only
## By group, tolerance slope as a function of methylation residuals:
modFULL <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*No.Worms + (1|brotherPairID) + (1|Sex),
                data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"),])
## Model selection:
step(modFULL, reduce.random = F) # Model found: BCI ~ (1 | brotherPairID) + (1 | Sex)

modFULL <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*PAT + (1|brotherPairID) + (1|Sex),
                data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"),])
## Model selection:
step(modFULL, reduce.random = F) # Model found: BCI ~ PAT + (1 | brotherPairID) + (1 | Sex)

## PCA based on all methylation values in G2
# heavy, run if needed
run = FALSE

if (run ==TRUE){
  # 1. get raw values
  percmeth = percMethylation(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP)
  
  # Run PCA on complete data (CpG covered in all fish)
  PCA_allpos <- myPCA(x = t(na.omit(percmeth)), incomplete = FALSE)
  
  # We perform a PCA on the CpG sites covered in all
  # G2 individuals. We first perform a test on the complete dataset.
  
  ### PCA plot with associated colors for treatments
  fviz_pca_ind(PCA_allpos$res.PCA,  label="none", habillage=PCA_allpos$metadata$trtG1G2) +
    scale_color_manual(values = colOffs)+
    scale_shape_manual(values=c(19,19,19,19))
  
  ### PCA plot with associated colors for brother pair
  fviz_pca_ind(PCA_allpos$res.PCA,  label="none", habillage=as.factor(PCA_allpos$metadata$brotherPairID))
  
  # The function dimdesc() can be used to identify the most correlated variables with a given principal component.
  mydimdesc <- dimdesc(PCA_allpos$res.PCA, axes = c(1,2), proba = 0.05)
  
  nrow(mydimdesc$Dim.1$quanti) #CpG sites most correlated (p < 0.05) with the first principal component , and `r nrow(mydimdesc$Dim.2$quanti)` with the second principal component.
  # The 2 first PCA axes do not explain BCI (p<0.05)
  
  ### How much of the BCI variance is explained by each variables?
  
  # Percentage of variance explained by each factor:
  formula(PCA_allpos$modSel) # BCI ~ No.Worms + PAT + (1 | brotherPairID) + (1 | Sex) + No.Worms:PAT
  mod_noworms = lmer(BCI ~ PAT + (1 | brotherPairID) + (1 | Sex), data = PCA_allpos$metadata)
  mod_noPAT = lmer(BCI ~ No.Worms + (1 | brotherPairID) + (1 | Sex), data = PCA_allpos$metadata)
  
  # R2c conditional R2 value associated with fixed effects plus the random effects.
  A = (MuMIn::r.squaredGLMM(PCA_allpos$modSel)[2] -
         MuMIn::r.squaredGLMM(mod_noworms)[2])*100
  
  B = (MuMIn::r.squaredGLMM(PCA_allpos$modSel)[2] -
         MuMIn::r.squaredGLMM(mod_noPAT)[2])*100
  
  round(A, 2)#% of the variance in associated with the parasite load (number of worms)
  round(B, 2)#% of the variance in associated with the paternal treatment
}
## Pretty summary picture
# Set up scatterplot
scatterplot <- ggplot(fullMetadata_OFFS,
                      aes(x = res_Nbr_methCpG_Nbr_coveredCpG,
                          y = BCI, fill=trtG1G2)) +
  geom_point(pch=21, size =3, alpha = .8) +
  guides(color = "none") +
  scale_fill_manual(values = colOffs, name = "Treatment",
                    labels = c("G1 control - G2 control", "G1 control - G2 exposed", "G1 exposed - G2 control", "G1 exposed - G2 exposed")) +
  theme(plot.margin = margin()) + theme_bw() +
  theme(legend.position = "none") +
  xlab("Methylation residuals (methylated sites/coverage")+
  ylab("Body Condition Index")

# Define marginal histogram
marginal_distribution <- function(x, var, group) {
  ggplot(x, aes_string(x = var, fill = group)) +
    # geom_histogram(bins = 30, alpha = 0.4, position = "identity") +
    geom_density(alpha = 0.6, size = 0.2) +
    guides(fill = "none") +
    scale_fill_manual(values = colOffs) +
    theme_void() +
    theme(plot.margin = margin())
}

# Set up marginal histograms
x_hist <- marginal_distribution(fullMetadata_OFFS, "res_Nbr_methCpG_Nbr_coveredCpG", "trtG1G2")
y_hist <- marginal_distribution(fullMetadata_OFFS, "BCI", "trtG1G2") +
  coord_flip()

# Align histograms with scatterplot
aligned_x_hist <- align_plots(x_hist, scatterplot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, scatterplot, align = "h")[[1]]

# Arrange plots
cowplot::plot_grid(
  aligned_x_hist, NULL, scatterplot, aligned_y_hist, ncol = 2, nrow = 2, rel_heights = c(0.2, 1), rel_widths = c(1, 0.2)
)
