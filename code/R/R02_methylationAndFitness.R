## Analyses
## A. Balard
## November 2021

library(plyr) # for join (keep row order)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)
library(ggsignif) ## for significance bars on ggplot
library(lme4) ## for mixed models
library(nlme) ## for mixed models
library(tidyverse)
library(emmeans) ## for post-hoc Tukey tests

## load custom functions
source("customRfunctions.R")

## Load samples metadata
source("R01.3_loadMetadata.R")

## Load previously united methylkit data

## define in which machine we're working (apocrita or mythinkpad)
machine="apocrita"
## machine="mythinkpad"
source("R01.4_loadMethyldata.R")
## Data getting loaded:
# uniteCovALL_woSexAndUnknowChr -> CpG positions shared by all fish
# uniteCovALL_G1_woSexAndUnknowChr-> CpG positions shared by all parents
# uniteCovALL_G2_woSexAndUnknowChr -> CpG positions shared by all offprings
# uniteCov6_G1_woSexAndUnknowChr -> CpG positions shared by half the parents in each trt group
# uniteCov14_G2_woSexAndUnknowChr -> CpG positions shared by half the offspring in each trt group

#####################################################################
## Compare fitness traits between the different offsprings groups ###
## Follow up of Sagonas 2020 & Ferre Ortega's master's dissertation #
#####################################################################

## Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).
fullMetadata_OFFS$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|Family), data=fullMetadata_OFFS))

## and for parents (no sex difference, only males):
fullMetadata_PAR$BCI <- residuals(lmer(Wnettofin ~ Slfin + (1|Family), data=fullMetadata_PAR))

## Effect of paternal treatment on body condition of offspring:
## Kaufmann et al. 2014:
# To investigate in which way paternal G1 exposure affected offspring tolerance,
# we tested how the relationship between G2 body condition and infection intensity
# was affected by paternal G1 exposure. This was tested in a linear mixed model on 
# G2 body condition with paternal G1 treatment and the interaction between 
# paternal G1 treatment and G2 infection intensity as fixed effects. Maternal 
# half-sibship identity was set as a random effect

## Effect of paternal exposure on tolerance:
modTol <- lme(BCI ~ patTrt + patTrt:No.Worms,random=~1|Family,data=fullMetadata_OFFS)
anova(modTol)

## Or modTol <- lme(BCI ~ patTrt*No.Worms,random=~1|Family,data=fullMetadata_offs)

myBCdf <- fullMetadata_OFFS %>% group_by(patTrt, No.Worms) %>% 
  summarise(BCI = mean(BCI)) %>% data.frame()
ggplot(fullMetadata_OFFS, aes(x=No.Worms, y = BCI, group = patTrt, col = patTrt))+
  geom_point() + geom_line(data=myBCdf)+
  geom_point(data=myBCdf, aes(fill = patTrt), col = "black", size = 3, pch = 21)+
  scale_color_manual(values = c("gray", "red"))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_bw()

## Effect of treatment groups of offspring on body condition:
## Kaufmann et al. 2014:
# The linear mixed effect model (nlme function in R) included G2 body condition as dependent variable, 
# sex, G2 treatment (exposed vs. control), paternal G1 treatment (exposed vs. control) 
# and their interactions as fixed effects as well as maternal G2 half-sibship identity as a random effect
mod1 <- lme(BCI ~ offsTrt * patTrt, random=~1|Family,data=fullMetadata_OFFS)
anova(mod1) # strong significant effect of both offspring trt & paternal + interactions

mod1.2 <- lme(BCI ~  trtG1G2, random=~1|Family,data=fullMetadata_OFFS)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than 
## every other group, same as Kaufmann et al. 2014

ggplot(fullMetadata_OFFS, aes(x=trtG1G2, y = BCI))+
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
  theme_bw()

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
## Calculate (1) Mfr per site, (2) biallelic equivalent, (3) epi-Fst and epi-FIS and (4) number of methylated sites in uniteCov2 (CpG shared by all fish, after filtering and normalising)

## Part 1: check correlation between methylation presence (beta > 0.7) and coverage in our 5 datasets
mycheckCorFUN <- function(myUniteCov, myMetaData){
    ## (1) Mfr per site
    df <- methylKit::getData(myUniteCov)
    MfrData <- data.frame(matrix(ncol=nrow(myMetaData),
                                 nrow=nrow(df)))
    namevector <- paste0("Mfr", 1:nrow(myMetaData))
    vector <- 1:nrow(myMetaData)
    for(i in vector){
        colnames(MfrData)[i] <- namevector[i]
        MfrData[i] <- df[paste0("numCs", i)]/
            df[paste0("coverage", i)]
    }
    ## (2) biallelic equivalent
### Sagonas et al. 2020 MBE " A  number  of  methylated  sites/regions  were  esti-mated by converting the MFr into ordinal data: sites/regionswith little or no methylation (MFr<30%) were annotated as0  and  treated  as  no  methylated  sites/regions,  sites/regionswith  intermediate  methylation  levels  (30%<MFr<70%)were considered as heterozygote sites/regions and convertedinto  1,  whereas  sites/regions  with  high  or  fixed  methylation(MFr>70%) were treated as homozygous at this site/regionsand were annotated as 2."

    getBiallVal <- function(x){
        y=NA
        if(x <= 0.3 & !is.na(x)){
            y = 0
        } else if(x > 0.3 & x < 0.7 & !is.na(x)){
            y = 1
        } else if(x >= 0.7 & !is.na(x)){
            y = 2
        }
        return(y)
    }
    getBiallValVec <- function(vec){sapply(vec, getBiallVal)}
    system.time(
        MbiallData <- sapply(MfrData, function(x) {getBiallValVec(x)})
    )

    ## rename correctly the full biallelic dataframe:
    colnames(MbiallData) <- myUniteCov@sample.ids

    ## make a data.frame with the sample ID and the biallelic methylation inf
    ## To target methylated regions, we selected only sites in which methylation  ratio  was  over  70%
    dfBiallMeth <- data.frame(
        biallmethy = apply(MbiallData, 2, function(v) sum(v == 2, na.rm=TRUE)),
        ID=myUniteCov@sample.ids)

    dfBiallMeth

    ## Append myMetaData with this new info
    idorder <- order(as.numeric(gsub("S", "", myMetaData$ID)))
    myMetaData <- merge(myMetaData, dfBiallMeth, by = "ID")
    myMetaData <- myMetaData[idorder, ]
    
    return(myMetaData)
}

fullMetadata_ALL <- mycheckCorFUN(uniteCovALL_woSexAndUnknowChr, fullMetadata)

fullMetadata_PAR_ALL  <- mycheckCorFUN(uniteCovALL_G1_woSexAndUnknowChr, fullMetadata_PAR)

fullMetadata_PAR_half <- mycheckCorFUN(uniteCov6_G1_woSexAndUnknowChr, fullMetadata_PAR)

fullMetadata_OFFS_ALL  <- mycheckCorFUN(uniteCovALL_G2_woSexAndUnknowChr, fullMetadata_OFFS)

fullMetadata_OFFS_half  <- mycheckCorFUN(uniteCov14_G2_woSexAndUnknowChr, fullMetadata_OFFS)

############
## Question: in our 5 datasets,
## Is there a correlation between nbr of methylated sites and coverage?
cor.test(fullMetadata_ALL$biallmethy,fullMetadata_ALL$M.Seqs_rawReads, method="pearson")
##t = 0.94062, df = 133, p-value = 0.3486

cor.test(fullMetadata_PAR_ALL$biallmethy,
         fullMetadata_PAR_ALL$M.Seqs_rawReads, method="pearson")
##t = -0.085834, df = 22, p-value = 0.9324

cor.test(fullMetadata_PAR_half$biallmethy,
         fullMetadata_PAR_half$M.Seqs_rawReads, method="pearson")
## t = 4.106, df = 22, p-value = 0.0004657 cor=0.6586753

cor.test(fullMetadata_OFFS_ALL$biallmethy,
         fullMetadata_OFFS_ALL$M.Seqs_rawReads, method="pearson")
## t = 1.0553, df = 109, p-value = 0.2936

cor.test(fullMetadata_OFFS_half$biallmethy,
         fullMetadata_OFFS_half$M.Seqs_rawReads, method="pearson")
##t = 4.1194, df = 109, p-value = 7.421e-05 cor=0.3670324

###########################################################
## Does RMS (ratio of methylated sites) changes with BCI ##
###########################################################
calcRMS <- function(myMetaData){
    myMetaData$biallmethy / (myMetaData$M.Seqs_rawReads*10e6)
}

fullMetadata_ALL$RMS <- calcRMS(fullMetadata_ALL)
fullMetadata_PAR_ALL$RMS <- calcRMS(fullMetadata_PAR_ALL)
fullMetadata_PAR_half$RMS <- calcRMS(fullMetadata_PAR_half)
fullMetadata_OFFS_ALL$RMS <- calcRMS(fullMetadata_OFFS_ALL)
fullMetadata_OFFS_half$RMS <- calcRMS(fullMetadata_OFFS_half)

## with family as random factor: not significant
## in parents:
anova(lme(RMS ~ BCI, random=~ 1|Family, data = fullMetadata_PAR_half))

## in offsprings:
anova(lme(RMS ~ BCI, random=~ 1|Family, data = fullMetadata_OFFS_half))
anova(lme(RMS ~ BCI*patTrt, random=~ 1|Family, data = fullMetadata_OFFS_half))

##############################################################################
## Does RMS (ratio of methylated sites) differ per trt group in offsprings? ##
##############################################################################

# simple lm, compare to null model: not signif
anova(lm(RMS ~ outcome, data = fullMetadata_OFFS_half)) 

# with family as random factor: not significant
anova(lme(RMS ~ outcome, random= ~1|Family, data = fullMetadata_OFFS_half))

# ggplot(fullMetadata_offs, aes(x = outcome, y = RMS, fill =outcome)) +
#     geom_boxplot() +
#     scale_fill_manual(values=c("grey","red")) +
#     theme_minimal()

## And per group? Not significant
mod <- lme(RMS ~ trtG1G2, random= ~1|Family, data = fullMetadata_OFFS)
anova(mod)

# ggplot(fullMetadata_offs, aes(x = trtG1G2, y = RMS)) +
#   geom_violin() +
#   geom_boxplot(width=.3) +
#   theme_minimal()

#########################################################################################
## Does RMS (ratio of methylated sites) correlated with parasite load in infected fish? ##
## -> NOPE
modRMS <- lme(RMS ~ patTrt*No.Worms, random=~1|Family, data=fullMetadata_OFFS_half[fullMetadata_OFFS_half$No.Worms > 0,])
anova(modRMS)

#plot(modRMS)

myRMSdf <- fullMetadata_OFFS_half[fullMetadata_OFFS_half$No.Worms > 0,] %>% group_by(patTrt, No.Worms) %>% 
  summarise(RMS = mean(RMS)) %>% data.frame()
myRMSdf

# ggplot(fullMetadata_offs[fullMetadata_offs$No.Worms > 0,], aes(x=No.Worms, y = RMS, group = patTrt, col = patTrt))+
#   geom_point() + geom_line(data=myRMSdf)+
#   geom_point(data=myRMSdf, aes(fill = patTrt), col = "black", size = 3, pch = 21)+
#   scale_color_manual(values = c("gray", "red"))+
#   scale_fill_manual(values = c("gray", "red"))+
#   theme_bw()

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
