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

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
## For further analyses: CpG covered in at least 2 then 6 individuals per group
load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov2_woSexAndUnknownChr.RData")
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCov6_woSexAndUnknowChr.RData")

## Load samples metadata
## Metylation metadata merged with Kostas files
fullMetadata <- read.csv("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/fullMetadata137_Alice.csv")
## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))
## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)

#####################################################################
## Compare fitness traits between the different offsprings groups ###
## Follow up of Sagonas 2020 & Ferre Ortega's master's dissertation #
#####################################################################
fullMetadata_offs <- fullMetadata[fullMetadata$Generat %in% "O",]
fullMetadata_offs$trtG1G2 <- droplevels(fullMetadata_offs$trtG1G2)

## Create variable for offsping and parents separated
fullMetadata_offs$offsTrt <- "controlO"
fullMetadata_offs$offsTrt[fullMetadata_offs$Tr %in% c("TT", "CT")] <- "infectedO"
fullMetadata_offs$patTrt <- "controlP"
fullMetadata_offs$patTrt[fullMetadata_offs$Tr %in% c("TC", "TT")] <- "infectedP"
## Sanity check
table(fullMetadata_offs$offsTrt, fullMetadata_offs$trtG1G2)
table(fullMetadata_offs$patTrt, fullMetadata_offs$trtG1G2)

## Kaufmann et al. 2014: Body condition of the G2 fish, an estimate of fish health and a predictor
# of energy reserves and reproductive success, was calculated using there residuals from the 
# regression of body mass on body length (Chellappaet al.1995).
fullMetadata_offs$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|Family), data=fullMetadata_offs))

## Effect of paternal treatment on body condition of offspring:
## Kaufmann et al. 2014:
# To investigate in which way paternal G1 exposure affected offspring tolerance,
# we tested how the relationship between G2 body condition and infection intensity
# was affected by paternal G1 exposure. This was tested in a linear mixed model on 
# G2 body condition with paternal G1 treatment and the interaction between 
# paternal G1 treatment and G2 infection intensity as fixed effects. Maternal 
# half-sibship identity was set as a random effect

## Effect of paternal exposure on tolerance:
modTol <- lme(BCI ~ patTrt + patTrt:No.Worms,random=~1|Family,data=fullMetadata_offs)
anova(modTol)

## Or modTol <- lme(BCI ~ patTrt*No.Worms,random=~1|Family,data=fullMetadata_offs)

myBCdf <- fullMetadata_offs %>% group_by(patTrt, No.Worms) %>% 
  summarise(BCI = mean(BCI)) %>% data.frame()
ggplot(fullMetadata_offs, aes(x=No.Worms, y = BCI, group = patTrt, col = patTrt))+
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
mod1 <- lme(BCI ~ offsTrt * patTrt, random=~1|Family,data=fullMetadata_offs)
anova(mod1) # strong significant effect of both offspring trt & paternal + interactions

mod1.2 <- lme(BCI ~  trtG1G2, random=~1|Family,data=fullMetadata_offs)
## pairwise posthoc test
emmeans(mod1.2, list(pairwise ~ trtG1G2), adjust = "tukey")
## Control father - treatment offspring has a strongly significantly lower BC than 
## every other group, same as Kaufmann et al. 2014

ggplot(fullMetadata_offs, aes(x=trtG1G2, y = BCI))+
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

## Nbr samples: 137
nrow(fullMetadata)
# Mean nbr of million reads: 11.2
mean(fullMetadata$M.Seqs_rawReads)
# 95% confidence interval: 0.35
qnorm(0.975)*sd(fullMetadata$M.Seqs_rawReads)/sqrt(nrow(fullMetadata))
# Average mapping efficiency +/-SD = 85.5% +/-0.47
mean(fullMetadata$MappingEfficiency.BSBoldvsGynogen)
qnorm(0.975)*sd(fullMetadata$MappingEfficiency.BSBoldvsGynogen)/sqrt(nrow(fullMetadata))

########################
## Calculate (1) Mfr per site, (2) biallelic equivalent, (3) epi-Fst and epi-FIS and (4) number of methylated sites in uniteCov2 (CpG shared by at least 2 fish, after filtering and normalising)

## (1) Mfr per site
reRun= FALSE
if(reRun == TRUE){
    df <- methylKit::getData(uniteCov2_woSexAndUnknowChr)
    MfrData2 <- data.frame(matrix(ncol=137, nrow=nrow(df)))
    namevector <- paste0("Mfr", 1:137)
    vector <- 1:137
    for(i in vector){
        colnames(MfrData2)[i] <- namevector[i]
        MfrData2[i] <- df[paste0("numCs", i)]/df[paste0("coverage", i)]
    }
}

## (2) biallelic equivalent
### Sagonas et al. 2020 MBE " A  number  of  methylated  sites/regions  were  esti-mated by converting the MFr into ordinal data: sites/regionswith little or no methylation (MFr<30%) were annotated as0  and  treated  as  no  methylated  sites/regions,  sites/regionswith  intermediate  methylation  levels  (30%<MFr<70%)were considered as heterozygote sites/regions and convertedinto  1,  whereas  sites/regions  with  high  or  fixed  methylation(MFr>70%) were treated as homozygous at this site/regionsand were annotated as 2."

rerun= FALSE
if(rerun == TRUE){
    print("long part started, needs ~10min")
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
        MbiallData2 <- sapply(MfrData2, function(x) {getBiallValVec(x)})
    )
#### Saving point ####
    save(MfrData2, MbiallData2, file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/MethylationFrequency_2fishmin.RData")
    print("long truc saved")
    stop("We stop here for now")
    print("Script did NOT end!")
}

#### Load data : MfrData2 & MbiallData2 ####
#load(file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/MethylationFrequency_2fishmin.RData")

## rename correctly the full biallelic dataframe:
colnames(MbiallData2) <- uniteCov2_woSexAndUnknowChr@sample.ids

## make a data.frame with the sample ID and the biallelic methylation inf
## To target methylated regions, we selectedonly sites in which methylation  ratio  was  over  70%
dfBiallMeth <- data.frame(
    biallmethy = apply(MbiallData2, 2, function(v) sum(v == 2, na.rm=TRUE)),
    ID=uniteCov2_woSexAndUnknowChr@sample.ids)

dfBiallMeth

## Append fullMetadata with this new info
fullMetadata <- merge(fullMetadata, dfBiallMeth, sort=F)
fullMetadata_offs <- merge(fullMetadata_offs, dfBiallMeth, sort=F) # N = 113 offsp.

# Is there a correlation between nbr of methylated sites and coverage? YES
cor.test(fullMetadata$biallmethy,fullMetadata$M.Seqs_rawReads, method="pearson")
## t = 6.8881, df = 135, p-value = 1.954e-10 cor 0.5099544

###########################################################
## Does RMS (ratio of methylated sites) changes with BCI ##
###########################################################
fullMetadata$RMS <- fullMetadata$biallmethy / (fullMetadata$M.Seqs_rawReads*10e6)
fullMetadata_offs$RMS <- fullMetadata_offs$biallmethy / (fullMetadata_offs$M.Seqs_rawReads*10e6)

## with family as random factor: not significant
mod <- lme(RMS ~ BCI, random=~ 1|Family, data = fullMetadata_offs)
anova(mod)

mod <- lme(RMS ~ BCI*patTrt, random=~ 1|Family, data = fullMetadata_offs)
anova(mod)

##############################################################################
## Does RMS (ratio of methylated sites) differ per trt group in offsprings? ##
##############################################################################


# simple lm, compare to null model: not signif
anova(lm(RMS ~ outcome, data = fullMetadata_offs)) 

# with family as random factor: not significant
mod <- lme(RMS ~ outcome, random= ~1|Family, data = fullMetadata_offs)
anova(mod)

# ggplot(fullMetadata_offs, aes(x = outcome, y = RMS, fill =outcome)) +
#     geom_boxplot() +
#     scale_fill_manual(values=c("grey","red")) +
#     theme_minimal()

## And per group? Not significant
mod <- lme(RMS ~ trtG1G2, random= ~1|Family, data = fullMetadata_offs)
anova(mod)

# ggplot(fullMetadata_offs, aes(x = trtG1G2, y = RMS)) +
#   geom_violin() +
#   geom_boxplot(width=.3) +
#   theme_minimal()

#########################################################################################
## Does RMS (ratio of methylated sites) correlated with parasite load in infected fish? ##
## -> NOPE
modRMS <- lme(RMS ~ patTrt*No.Worms, random=~1|Family, data=fullMetadata_offs[fullMetadata_offs$No.Worms > 0,])
anova(modRMS)

#plot(modRMS)

myRMSdf <- fullMetadata_offs[fullMetadata_offs$No.Worms > 0,] %>% group_by(patTrt, No.Worms) %>% 
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
library(genepop)

library(adegenet)

test <- head(MbiallData2, 15)

offspringID <- fullMetadata_offs$ID

offspring_gen = df2genind(test[,colnames(test) %in% offspringID],
                          ploidy = 2, ind.names = offspringID,
                          pop = fullMetadata_offs$trtG1G2, sep = "")


seafan_gen = import2genind("Pinkseafan_13MicrosatLoci.gen", ncode = 3, quiet = TRUE)

head(seafan_gen)



################# PCA
## PCA analysis on our samples: plot a scree plot for importance of components
p=PCASamples(uniteCovALL_N137_final, screeplot=TRUE, obj.return = T) # first axis very important
s=summary(p)
#create scree plot
library(ggplot2)
qplot(c(1:137), s$importance[2,]) + 
  geom_line() + 
  xlab("Principal Component") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L), name = "Variance Explained")+
  ggtitle("Scree Plot") +
  theme_bw()

## PCA: validate the absence of strong sex influence on the methylation
uniteCovALL_N137_final_SEX = uniteCovALL_N137_final
uniteCovALL_N137_final_SEX@treatment = as.numeric(as.factor(metadata$Sex))
PCASamples(uniteCovALL_N137_final_SEX); title(sub="colored by sex")
PCASamples(uniteCovALL_N137_final_SEX, comp = c(3,4)); title(sub="colored by sex")

## PCA: check the Family influence on the methylation pattern
uniteCovALL_N137_final_FAM = uniteCovALL_N137_final
uniteCovALL_N137_final_FAM@treatment = as.numeric(as.factor(metadata$Family))
PCASamples(uniteCovALL_N137_final_FAM); title(sub="colored by family")
PCASamples(uniteCovALL_N137_final_FAM, comp = c(3,4)); title(sub="colored by family")

## PCA: check the Pattern treatment influence on the methylation pattern
uniteCovALL_N137_final_PAT = uniteCovALL_N137_final
metadata$PAT="Exposed father group"
metadata$PAT[metadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"
uniteCovALL_N137_final_PAT@treatment = as.numeric(as.factor(metadata$PAT))
PCASamples(uniteCovALL_N137_final_PAT); title(sub="colored by paternal exposure group")
PCASamples(uniteCovALL_N137_final_PAT, comp=c(3,4)); title(sub="colored by paternal exposure group")



