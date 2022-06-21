## Analyses
## A. Balard
## November 2021

machine="mythinkpad" # define the machine we work on
loadALL = FALSE # only load CpG shared by half fish per trt group
source("R02.3_DATALOAD.R")

## Data getting loaded:
uniteCovALL_woSexAndUnknowChr #-> 55530 CpG positions shared by all fish
uniteCov6_G1_woSexAndUnknowChrOVERLAP #-> 1001880 CpG positions shared by half the parents in each trt group, overlapping with the parental ones
uniteCov14_G2_woSexAndUnknowChrOVERLAP#-> 1001880 CpG positions shared by half the offspring in each trt group, overlapping with the offspring ones

#####################################################################
## Compare fitness traits between the different offsprings groups ###
## Follow up of Sagonas 2020 & Ferre Ortega's master's dissertation #
#####################################################################

## Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).
fullMetadata_OFFS$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|brotherPairID), data=fullMetadata_OFFS))

## and for parents (no sex difference, only males):
fullMetadata_PAR$BCI <- residuals(lmer(Wnettofin ~ Slfin + (1|brotherPairID), data=fullMetadata_PAR))

## Effect of paternal treatment on body condition of offspring:
## Kaufmann et al. 2014:
# To investigate in which way paternal G1 exposure affected offspring tolerance,
# we tested how the relationship between G2 body condition and infection intensity
# was affected by paternal G1 exposure. This was tested in a linear mixed model on 
# G2 body condition with paternal G1 treatment and the interaction between 
# paternal G1 treatment and G2 infection intensity as fixed effects. Maternal 
# half-sibship identity was set as a random effect

############################################
## Effect of paternal exposure on tolerance:
modTol <- lmer(BCI ~ patTrt + patTrt:No.Worms + (1|brotherPairID), data=fullMetadata_OFFS)

stepcAIC(modTol)

pred <- ggpredict(modTol, terms = c("No.Worms", "patTrt"))
plot(pred, add.data = TRUE)
## The slope of BCI on nbrworms varies upon parental treatment = parental treatment varies with tolerance

## Effect of treatment groups of offspring on body condition:
## Kaufmann et al. 2014:
# The linear mixed effect model (nlme function in R) included G2 body condition as dependent variable, 
# sex, G2 treatment (exposed vs. control), paternal G1 treatment (exposed vs. control) 
# and their interactions as fixed effects as well as maternal G2 half-sibship identity as a random effect
mod1 <- lme(BCI ~ offsTrt * patTrt, random=~1|brotherPairID,data=fullMetadata_OFFS)
anova(mod1) # strong significant effect of both offspring trt & paternal + interactions

mod1.2 <- lme(BCI ~  trtG1G2, random=~1|brotherPairID,data=fullMetadata_OFFS)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than 
## every other group, same as Kaufmann et al. 2014

ggplot(fullMetadata_OFFS, aes(x=trtG1G2, y = BCI, fill=trtG1G2))+
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
  theme_bw() + theme(legend.position = "none")

############################################
## Calculate tolerance as a reaction norm ##
############################################
mod1 <- lmer(BCI ~ No.Worms : trtG1G2 + (1|brotherPairID) + (1|Sex), 
             data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"),], 
             REML = FALSE)
mod0 <- lmer(BCI ~ No.Worms + (1|brotherPairID) + (1|Sex),
             data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"),], 
             REML = FALSE)

lrtest(mod0, mod1)
# LRT: parental treatment groups (infected or not): G = 18.598, df = 6, p = 1.614e-05 ***

modFULL <- lmer(BCI ~ No.Worms : trtG1G2 + (1|brotherPairID) + (1|Sex), 
                data = fullMetadata_OFFS[fullMetadata_OFFS$trtG1G2 %in% c("NE_exposed", "E_exposed"),])

coef(modFULL)
# infected father: BCI = 7.103 NoWorms + various intercepts by family and sex
# control father: BCI = -15.594 NoWorms + various intercepts by family and sex

# plot fixed effects
pred <- ggpredict(modFULL, terms = c("No.Worms", "trtG1G2"))

plot(pred, add.data = TRUE)+ theme_bw() +
  scale_x_continuous(name = "Number of parasites", breaks = 1:10)+
  scale_y_continuous("Body Condition Index")+
  ggtitle("Predicted values of Body Condition Index")+
  scale_color_manual(values = colOffs[c(2,4)])

## And in fathers:
modFULL <- lmer(BCI ~ No.Worms + (1|Family), 
                data = fullMetadata_PAR[fullMetadata_PAR$trtG1G2 %in% "Exposed",])

coef(modFULL)
# BCI = -7.3 NoWorms + 38.5

# plot fixed effects
pred <- ggpredict(modFULL, terms = "No.Worms")

plot(pred, add.data = TRUE)+ theme_bw() +
  scale_x_continuous(name = "Number of parasites", breaks = 1:10)+
  scale_y_continuous("Body Condition Index")+
  ggtitle("Predicted values of Body Condition Index")

#######################################################
## Nbr/Ratio of Methylated Sites in different groups ##
#######################################################

## Nbr samples: 135
nrow(fullMetadata)
# Mean nbr of million reads: 11.3
mean(fullMetadata$M.Seqs_rawReads)
# 95% confidence interval: 0.33
qnorm(0.975)*sd(fullMetadata$M.Seqs_rawReads)/sqrt(nrow(fullMetadata))
# Average mapping efficiency +/-SD = 85.4% +/-0.48
mean(fullMetadata$MappingEfficiency.BSBoldvsGynogen)
qnorm(0.975)*sd(fullMetadata$MappingEfficiency.BSBoldvsGynogen)/sqrt(nrow(fullMetadata))

########################
## Calculate (1) number of methylated sites and 
## (2) ratio of methylated sites (to account for coverage bias)

mycalcRMS <- function(myUniteCov, myMetaData){
  percMethMat = methylKit::percMethylation(myUniteCov)
  # create a dataframe with all info
  percMethDF = data.frame(SampleID = colnames(percMethMat),
                          ## number of methylated sites
                          Nbr_methCpG = colSums(percMethMat>=70 & !is.na(percMethMat)),#512493
                          ## number of sites covered in this sample
                          Nbr_coveredCpG = colSums(!is.na(percMethMat)), #1015735
                          ## number of sites NOT covered in this sample
                          Nbr_NOTcoveredCpG = colSums(is.na(percMethMat)))
  ## RMS in this sample based on covered sites
  percMethDF$RMS_coveredCpG = percMethDF$Nbr_methCpG / percMethDF$Nbr_coveredCpG #0.5045538
  ## merge with original metadata:
  myMetaData = merge(myMetaData, percMethDF)
  # calculate also RMS global, considering CpG covered or not (to compare)
  myMetaData$RMS_allCpG_coveredOrNot = myMetaData$Nbr_methCpG / (myMetaData$M.Seqs_rawReads*10e6)
  # calculate residuals of nbr of methCpG by nbr of covered CpG
  myMetaData$res_Nbr_methCpG_Nbr_coveredCpG = residuals(
    lm(myMetaData$Nbr_methCpG ~ myMetaData$Nbr_coveredCpG))
  return(myMetaData)
}

fullMetadata_ALL <- mycalcRMS(uniteCovALL_woSexAndUnknowChr, fullMetadata)

fullMetadata_PAR_half <- mycalcRMS(uniteCov6_G1_woSexAndUnknowChrOVERLAP, fullMetadata_PAR)

fullMetadata_OFFS_half  <- mycalcRMS(uniteCov14_G2_woSexAndUnknowChrOVERLAP, fullMetadata_OFFS)

## Question: in our 5 datasets,
## Is there a correlation between nbr of methylated sites and coverage? 
cor.test(fullMetadata_PAR_half$Nbr_coveredCpG,
         fullMetadata_PAR_half$Nbr_methCpG, method = "spearman")
## S = 350, p-value = 2.15e-06, rho = 0.85
ggplot(fullMetadata_PAR_half, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

## Check after RMS correction for coverage bias: CORRECTED (p-value = 0.4485)
cor.test(fullMetadata_PAR_half$Nbr_coveredCpG,
         fullMetadata_PAR_half$RMS_coveredCpG, method = "spearman")
ggplot(fullMetadata_PAR_half, aes(x=Nbr_coveredCpG, y=RMS_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

## and with residuals: COMPLETELY CORRECTED p-value = 0.9562
cor.test(fullMetadata_PAR_half$Nbr_coveredCpG,
         fullMetadata_PAR_half$res_Nbr_methCpG_Nbr_coveredCpG, method = "spearman")
ggplot(fullMetadata_PAR_half, aes(x=Nbr_coveredCpG, y=res_Nbr_methCpG_Nbr_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "Parents, CpG shared by half fish/trt")

############
## Offspring:
cor.test(fullMetadata_OFFS_half$Nbr_coveredCpG,
         fullMetadata_OFFS_half$Nbr_methCpG, method = "spearman")
## S = 20254, p-value < 2.2e-16 rho = 0.91
ggplot(fullMetadata_OFFS_half, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  scale_x_continuous("Number of cytosines covered") +
  scale_y_continuous("Number of methylated cytosines") +
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## Plot distance to residuals:
fit <- lm(Nbr_methCpG ~ Nbr_coveredCpG, data = fullMetadata_OFFS_half)
plotdf <- fullMetadata_OFFS_half
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
cor.test(fullMetadata_OFFS_half$Nbr_coveredCpG,
         fullMetadata_OFFS_half$RMS_coveredCpG, method = "spearman")
ggplot(fullMetadata_OFFS_half, aes(x=Nbr_coveredCpG, y=RMS_coveredCpG))+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  geom_smooth(method = "lm", col="black")+
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

## and with residuals: COMPLETELY CORRECTED p-value = 0.51
cor.test(fullMetadata_OFFS_half$Nbr_coveredCpG,
         fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG, method = "spearman")
ggplot(fullMetadata_OFFS_half, aes(x=Nbr_coveredCpG, y=res_Nbr_methCpG_Nbr_coveredCpG))+
  geom_point(aes(col=trtG1G2), size = 3)+ scale_color_manual(values = colOffs) +
  geom_smooth(method = "lm", col="black")+
  scale_x_continuous("Number of cytosines covered") +
  scale_y_continuous("Residuals of number of methylated cytosines\n on number of cytosines covered") +
  theme_bw() + ggtitle(label = "Offspring, CpG shared by half fish/trt")

################################
## Should we correct for sex? ##
################################

################
## Does Sex affect the number of methylated sites? YES
## + family as random factor
modFull <- lmer(Nbr_methCpG ~ trtG1G2 * Sex + (1|brotherPairID), 
                data = fullMetadata_OFFS_half, REML = F) # REML =F for model comparison
mod_noSex <- lmer(Nbr_methCpG ~ trtG1G2 + (1|brotherPairID), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noTrt <- lmer(Nbr_methCpG ~ Sex + (1|brotherPairID), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noInteractions <- lmer(Nbr_methCpG ~ trtG1G2 + Sex + (1|brotherPairID), 
                           data = fullMetadata_OFFS_half, REML = F)

lrtest(modFull, mod_noSex) # sex is VERY VERY significant p = 0.001776 **
lrtest(modFull, mod_noTrt) # trt is signif p = 0.0208 *
lrtest(modFull, mod_noInteractions) # interactions are significant 0.0151 *

## Plot
ggplot(fullMetadata_OFFS_half, aes(trtG1G2, Nbr_methCpG, group=interaction(trtG1G2, Sex))) + 
  facet_grid(~Sex) +
  geom_violin() +
  geom_boxplot(aes(fill = trtG1G2), width = 0.2) +
  geom_jitter(width = .1, size = 1, pch = 21, fill = "white") + 
  scale_fill_manual(values = colOffs) +
  theme_bw()  + theme(legend.position = "none")

################
## Does Sex affect the residuals of nbr of methylated sites by nbr of sites? YES
## + family as random factor
modFull <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ trtG1G2 * Sex + (1|brotherPairID), 
                data = fullMetadata_OFFS_half, REML = F) # REML =F for model comparison
mod_noSex <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ trtG1G2 + (1|brotherPairID), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noTrt <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ Sex + (1|brotherPairID), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noInteractions <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ trtG1G2 + Sex + (1|brotherPairID), 
                           data = fullMetadata_OFFS_half, REML = F)

lrtest(modFull, mod_noSex) # sex is significant p = 0.0002124 ***
lrtest(modFull, mod_noTrt) # trt is not significant any longer 
lrtest(modFull, mod_noInteractions) # interactions are are not significant any longer

## Plot
ggplot(fullMetadata_OFFS_half, aes(trtG1G2, res_Nbr_methCpG_Nbr_coveredCpG,
                                   group=interaction(trtG1G2, Sex))) + 
  facet_grid(~Sex) +
  geom_violin() +
  geom_boxplot(aes(fill = trtG1G2), width = 0.2) +
  geom_jitter(width = .1, size = 1, pch = 21, fill = "white") + 
  scale_fill_manual(values = colOffs) +
  theme_bw() + theme(legend.position = "none")

########################################################################
## Are mean residuals meth sites different following tolerance slope? ##
########################################################################
fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG_div1000 <- (fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG)/1000

mod_Tol.Meth <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*No.Worms + (1|brotherPairID)+ (1|Sex), 
                     data=fullMetadata_OFFS_half, REML = F)
mod_Tol.Meth_Noint <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000+No.Worms + (1|brotherPairID) + (1|Sex),
                           data=fullMetadata_OFFS_half, REML = F)
lrtest(mod_Tol.Meth, mod_Tol.Meth_Noint)
## The slope of BCI on nbrworms varies upon parental treatment = methylation does NOT vary with tolerance

## And just among infected?
mod_Tol.Meth.inf <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000*No.Worms + (1|brotherPairID)+ (1|Sex), 
                     data=fullMetadata_OFFS_half[fullMetadata_OFFS_half$outcome %in% "infected",], REML = F)
mod_Tol.Meth_Noint.inf <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000+No.Worms + (1|brotherPairID) + (1|Sex),
                           data=fullMetadata_OFFS_half[fullMetadata_OFFS_half$outcome %in% "infected",], REML = F)
lrtest(mod_Tol.Meth.inf, mod_Tol.Meth_Noint.inf)
# also NOT in infected ofspring

pred <- ggpredict(mod_Tol.Meth, terms = c("No.Worms", "res_Nbr_methCpG_Nbr_coveredCpG_div1000"))
plot(pred, add.data = TRUE)

############################################################
## Are mean residuals meth sites different following BCI? ##
############################################################
mod_BCI.Meth <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000 + (1|brotherPairID)+ (1|Sex), 
                     data=fullMetadata_OFFS_half, REML = F)
mod_BCI.Meth.NULL <- lmer(BCI ~ 1 + (1|brotherPairID)+ (1|Sex), 
                     data=fullMetadata_OFFS_half, REML = F)
lrtest(mod_BCI.Meth, mod_BCI.Meth.NULL)
# NOT in all ofspring

## And just among infected?
mod_BCI.Meth_inf <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_div1000 + (1|brotherPairID)+ (1|Sex), 
                     data=fullMetadata_OFFS_half[fullMetadata_OFFS_half$outcome %in% "infected",], REML = F)
mod_BCI.Meth_inf.NULL <- lmer(BCI ~ 1 + (1|brotherPairID)+ (1|Sex), 
                          data=fullMetadata_OFFS_half[fullMetadata_OFFS_half$outcome %in% "infected",], REML = F)
lrtest(mod_BCI.Meth_inf, mod_BCI.Meth_inf.NULL)
# also NOT in infected ofspring





pred <- ggpredict(mod_Tol.Meth, terms = c("No.Worms", "res_Nbr_methCpG_Nbr_coveredCpG_div1000"))
plot(pred, add.data = TRUE)
## The slope of BCI on nbrworms varies upon parental treatment = methylation does NOT vary with tolerance



#####################################################################
## Are mean residuals meth sites different in NE_E and E_E groups? ##
#####################################################################

## By group, tolerance slope as a function of methylation residuals:
modFULL <- lmer(BCI ~ No.Worms : trtG1G2 + (1|brotherPairID) + (1|Sex), 
                data = fullMetadata_OFFS_half[fullMetadata_OFFS_half$trtG1G2 %in% c("NE_exposed", "E_exposed"),])

coef(modFULL)[1]
# infected father: BCI = 7.1 NoWorms + various intercepts by family and sex
# control father: BCI = -15.6 NoWorms + various intercepts by family and sex

predict <- ggpredict(modFULL, terms = c("No.Worms", "trtG1G2"))

slope_NE_exposed <- predict[predict$group %in% "NE_exposed",][2,2] - predict[predict$group %in% "NE_exposed",][1,2]
slope_E_exposed <- predict[predict$group %in% "E_exposed",][2,2] - predict[predict$group %in% "E_exposed",][1,2]

slope_NE_exposed

meanMeth_NE_exposed <- mean(fullMetadata_OFFS_half[fullMetadata_OFFS_half$trtG1G2 %in% c("NE_exposed"),"res_Nbr_methCpG_Nbr_coveredCpG"])
meanMeth_E_exposed <- mean(fullMetadata_OFFS_half[fullMetadata_OFFS_half$trtG1G2 %in% c("E_exposed"),"res_Nbr_methCpG_Nbr_coveredCpG"])

df <- data.frame(group = c("NE_exposed", "E_exposed"),
                 slope = c(slope_NE_exposed, slope_E_exposed),
                 meanMeth = c(meanMeth_NE_exposed, meanMeth_E_exposed))

ggplot(df, aes(x = slope, y = meanMeth, col =group))+
  geom_point() + theme_bw()

## Low tolerance = high methylation / high tolerance = low methylation

modTOLfamsex <- lm(BCI ~ No.Worms : trtG1G2 : brotherPairID : Sex, 
                   data = fullMetadata_OFFS_half[fullMetadata_OFFS_half$trtG1G2 %in% c("NE_exposed", "E_exposed"),])

coef(modTOLfamsex)[1]


predTOLfamsex <- ggpredict(modTOLfamsex, terms = c("trtG1G2", "brotherPairID", "Sex"))

plot(predTOLfamsex)

##################################################
## Do residuals meth sites differ per trt group ##
##################################################

#############START
mod1.2 <- lme(BCI ~  trtG1G2, random=~1|brotherPairID,data=fullMetadata_OFFS)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than 
## every other group, same as Kaufmann et al. 2014

ggplot(fullMetadata_OFFS, aes(x=trtG1G2, y = BCI, fill=trtG1G2))+
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
  theme_bw() + theme(legend.position = "none")
#############END


# NB: we put sex (in offspring) and family as random factors

## PARENTS
mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ trtG1G2 + (1|brotherPairID), 
             data = fullMetadata_PAR_half, REML = F)
mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ 1 + (1|brotherPairID), 
             data = fullMetadata_PAR_half, REML = F)
lrtest(mod1, mod0) # not significant in parents

## OFFSPRINGS
mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ trtG1G2 + (1|brotherPairID) + (1|Sex), 
             data = fullMetadata_OFFS_half, REML = F)
mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ 1 + (1|brotherPairID) + (1|Sex), 
             data = fullMetadata_OFFS_half, REML = F)
lrtest(mod1, mod0) # not significant in offspring

## Plot
ggplot(fullMetadata_OFFS_half, aes(x=trtG1G2, y = res_Nbr_methCpG_Nbr_coveredCpG, fill=trtG1G2))+
  geom_boxplot()+
  scale_x_discrete("Treatment")+
  scale_y_continuous("Residuals of number of methylated cytosines\n on number of cytosines covered") +
  scale_fill_manual(values = colOffs)+
  theme_bw() + theme(legend.position = "none")

##############################################################################
## Decompose: do residuals methylated sites change with PAR & OFF treatment ##
##############################################################################
## OFFSPRINGS
modFull <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ patTrt * outcome + (1|brotherPairID) + (1|Sex), 
                data = fullMetadata_OFFS_half, REML = F)
mod_noPAT <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ outcome + (1|brotherPairID) + (1|Sex), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noOFF <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ patTrt + (1|brotherPairID) + (1|Sex), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noInt <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ patTrt + outcome + (1|brotherPairID) + (1|Sex), 
                  data = fullMetadata_OFFS_half, REML = F)

lrtest(modFull, mod_noPAT) # not significant
lrtest(modFull, mod_noOFF) # not significant
lrtest(modFull, mod_noInt) # not significant

########################################################################################
## Decompose: do residuals methylated sites change with BCI & parasite load treatment ##
########################################################################################
## OFFSPRINGS
modFull <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI * No.Worms + (1|brotherPairID) + (1|Sex), 
                data = fullMetadata_OFFS_half, REML = F)
mod_noBCI <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ No.Worms + (1|brotherPairID) + (1|Sex), 
                  data = fullMetadata_OFFS_half, REML = F)
mod_noNo.Worms <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI + (1|brotherPairID) + (1|Sex), 
                       data = fullMetadata_OFFS_half, REML = F)
mod_noInt <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI + No.Worms + (1|brotherPairID) + (1|Sex), 
                  data = fullMetadata_OFFS_half, REML = F)

lrtest(modFull, mod_noBCI) # not significant
lrtest(modFull, mod_noNo.Worms) # not significant
lrtest(modFull, mod_noInt) # not significant

### JUST BCI:
# NB: we put sex (in offspring) and family as random factors

df = data.frame(BCI = c(fullMetadata_PAR_half$BCI, fullMetadata_OFFS_half$BCI),
                brotherPairID = c(fullMetadata_PAR_half$brotherPairID, fullMetadata_OFFS_half$brotherPairID),
                res_Nbr_methCpG_Nbr_coveredCpG = c(fullMetadata_PAR_half$res_Nbr_methCpG_Nbr_coveredCpG, fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG))

## ALL
mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI + (1|brotherPairID), 
             data = df, REML = F)
mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ 1 + (1|brotherPairID), 
             data = df, REML = F)
lrtest(mod1, mod0) # not significant

## PARENTS
mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI + (1|brotherPairID), 
             data = fullMetadata_PAR_half, REML = F)
mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ 1 + (1|brotherPairID), 
             data = fullMetadata_PAR_half, REML = F)
lrtest(mod1, mod0) # not significant in parents

## OFFSPRING
mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ BCI + (1|brotherPairID) + (1|Sex), 
             data = fullMetadata_OFFS_half, REML = F)
mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG ~ 1 + (1|brotherPairID)+ (1|Sex), 
             data = fullMetadata_OFFS_half, REML = F)
lrtest(mod1, mod0) # not significant in offspring

# plot fixed effects depending on random effects
pred <- ggpredict(mod1, terms = c("BCI", "brotherPairID", "Sex"), type = "random")
plot(pred, ci = F, add.data = TRUE)+
  ggtitle("Predicted values of global methylation")+
  scale_y_continuous("Residuals of number of methylated cytosines\n on number of cytosines covered") +
  scale_x_continuous("Body Condition Index")

# Créer le graphique
p + facet_grid(
  dose ~ supp, 
  labeller = labeller(dose = dose.labs, supp = supp.labs)
)


## Plotted in the other direction
plot(pred, ci = F, add.data = TRUE)+
  coord_flip()

## In the other direction:
### to avoid scaling error:
# fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG_DIVBY1000 <- fullMetadata_OFFS_half$res_Nbr_methCpG_Nbr_coveredCpG/1000
# 
# mod1 <- lmer(BCI ~ res_Nbr_methCpG_Nbr_coveredCpG_DIVBY1000 + (1|Family) + (1|Sex), 
#              data = fullMetadata_OFFS_half, REML = F)
# mod0 <- lmer(BCI ~ 1 + (1|Family)+ (1|Sex), 
#              data = fullMetadata_OFFS_half, REML = F)
# lrtest(mod1, mod0) # significant in offspring p = 0.0162 *
# 
# # plot fixed effects depending on random effects
# pred <- ggpredict(mod1, terms = c("res_Nbr_methCpG_Nbr_coveredCpG", "Family", "Sex"), type = "random")
# plot(pred, ci = F, add.data = TRUE)



##### Prettier picture:
# Set up scatterplot
scatterplot <- ggplot(fullMetadata_OFFS_half,
                      aes(x = res_Nbr_methCpG_Nbr_coveredCpG, 
                          y = BCI, fill=trtG1G2)) +
  geom_point(pch=21, size =3, alpha = .8) +
  guides(color = "none") +
  scale_fill_manual(values = colOffs) + 
  theme(plot.margin = margin()) +
  theme(legend.position = "none")

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
x_hist <- marginal_distribution(fullMetadata_OFFS_half, "res_Nbr_methCpG_Nbr_coveredCpG", "trtG1G2")
y_hist <- marginal_distribution(fullMetadata_OFFS_half, "BCI", "trtG1G2") +
  coord_flip()

# Align histograms with scatterplot
aligned_x_hist <- align_plots(x_hist, scatterplot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, scatterplot, align = "h")[[1]]

# Arrange plots
plot_grid(
  aligned_x_hist
  , NULL
  , scatterplot
  , aligned_y_hist
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)


#### TBC
#################################################################################
## Calcul of epi-FST/epi-FIS: Homogeneisation of methylation marks in infected? ##
################################################################################
# library(genepop)
# library(adegenet)
# 
# test <- head(MbiallData2, 15)
# 
# offspringID <- fullMetadata_OFFS$ID
# 
# offspring_gen = df2genind(test[,colnames(test) %in% offspringID],
#                           ploidy = 2, ind.names = offspringID,
#                           pop = fullMetadata_OFFS$trtG1G2, sep = "")
# 
# 
# seafan_gen = import2genind("Pinkseafan_13MicrosatLoci.gen", ncode = 3, quiet = TRUE)
# 
# head(seafan_gen)
