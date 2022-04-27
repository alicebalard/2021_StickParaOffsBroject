## Differential methylation analyses
## A. Balard
## February 2022

#################### Data load & preparation ####################
source("librariesLoading.R")
## load custom functions
source("customRfunctions.R")
## Load samples metadata
source("R02.1_loadMetadata.R")
## define in which machine we're working (apocrita or mythinkpad)
##machine="apocrita"
machine="mythinkpad"
## Load methylation data
loadALL = FALSE # only load CpG shared by half fish per trt group + ALL
source("R02.2_loadMethyldata.R")
## Load file containing length of each gynogen chromosomes
## grep "contig" gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff | awk '{print $1, $5}' > data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt
GYgynogff = read.table("../../data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff) = c("chrom","length")

###########################################
## Source the previously calculated DMS/DMR
## Parents (brotherPairID as covariates)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G1_half <- readRDS("../../data/DiffMeth/DMS15pc_G1_half.RDS"); nrow(DMS15pc_G1_half) # 3648
DMR15pc_G1_half <- readRDS("../../data/DiffMeth/DMR15pc_G1_half.RDS"); nrow(DMR15pc_G1_half) # 23
### DM from CpG positions shared by all the fish
DMS15pc_G1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G1_ALL.RDS"); nrow(DMS15pc_G1_ALL) # 125
# DMR15pc_G1_ALL returned 0 DMR

## Offspring (brotherPairID & Sex as covariates)
## Control G1 - G2(trt vs control)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_half.RDS")
nrow(DMS15pc_G2_controlG1_half) # 1197
DMR15pc_G2_controlG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_controlG1_half.RDS")
nrow(DMR15pc_G2_controlG1_half) # 6
### DM from CpG positions shared by all the fish
DMS15pc_G2_controlG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_controlG1_ALL.RDS")
nrow(DMS15pc_G2_controlG1_ALL) # 38
# DMR15pc_G2_controlG1_ALL returned 0 DMR

## Infected G1 - G2(trt vs control)
### DM from CpG positions shared by half the fish per trt
DMS15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_half.RDS")
nrow(DMS15pc_G2_infectedG1_half) # 690
DMR15pc_G2_infectedG1_half <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_half.RDS")
nrow(DMR15pc_G2_infectedG1_half) # 8
### DM from CpG positions shared by all the fish
DMS15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMS15pc_G2_infectedG1_ALL.RDS")
nrow(DMS15pc_G2_infectedG1_ALL) # 22
DMR15pc_G2_infectedG1_ALL <- readRDS("../../data/DiffMeth/DMR15pc_G2_infectedG1_ALL.RDS")
nrow(DMR15pc_G2_infectedG1_ALL) # 1
###########################################

###########################
## Function to get DMS info
myDMSinfo <- function(DMSobject, fromUniteCov){
  DMS = paste(DMSobject$chr, DMSobject$start, DMSobject$end)
  meth.diff = DMSobject$meth.diff
  direction = ifelse(DMSobject$meth.diff > 0, "hyper", "hypo")
  percentDMS = length(DMS)/nrow(fromUniteCov)*100
  return(list(DMS = DMS, meth.diff = meth.diff, direction = direction, percentDMS = percentDMS))
}

## Run the function
DMS_info_G1 <- myDMSinfo(DMS15pc_G1_half, uniteCov6_G1_woSexAndUnknowChrOVERLAP)
DMS_info_G2_G1c_final <- myDMSinfo(DMS15pc_G2_controlG1_half, uniteCov14_G2_woSexAndUnknowChrOVERLAP)
DMS_info_G2_G1i_final <- myDMSinfo(DMS15pc_G2_infectedG1_half,uniteCov14_G2_woSexAndUnknowChrOVERLAP)

## NB Kostas' results: "We found a total of 1,973 CpG sites out of 1,172,887 CpGs (0.17%)
# across the genome that showed at least 15% differential fractional methylation
# (differentially methylated site [DMS]; q < 0.01) between infected and uninfected fish"

## Here: number of CpG sites
nrow(uniteCov14_G2_woSexAndUnknowChrOVERLAP) # 1,001,880

## Parents comparison:
length(DMS_info_G1$DMS)# 3648 DMS
DMS_info_G1$percentDMS # 0.36% of the CpGs are DMS

## Offspring from control parents comparison:
length(DMS_info_G2_G1c_final$DMS) # 1197 DMS
DMS_info_G2_G1c_final$percentDMS # 0.12% of the CpGs are DMS

## Offspring from infected parents comparison:
length(DMS_info_G2_G1i_final$DMS) # 690 DMS
DMS_info_G2_G1i_final$percentDMS # 0.07% of the CpGs are DMS

####################################################################################
## Question: how are the beta values in the different groups for the parental DMS?##
####################################################################################
##############
## Prepare dataset
##############
PM_G1 <- getPMdataset(uniteCov = uniteCov6_G1_woSexAndUnknowChrOVERLAP, MD = fullMetadata_PAR, gener="parents")
PM_G2 <- getPMdataset(uniteCov = uniteCov14_G2_woSexAndUnknowChrOVERLAP, MD = fullMetadata_OFFS, gener="offspring")

head(PM_G1)
head(PM_G2)

table(fullMetadata_OFFS$trtG1G2, fullMetadata_OFFS$clutch.ID)

## What is the relative contribution of methylation to brother pair & paternal treatment?
## Test of VCA: variance component analysis https://cran.r-project.org/web/packages/VCA/vignettes/VCA_package_vignette.html
library(VCA)
myfitVCA <- fitVCA(form = BetaValue~(G1_trt+G2_trt)/brotherPairID/CpGSite,
                   Data = sample_n(PM_G2, 100))
print(myfitVCA, digits=4)

##############
## In parents
##############
parmod <- lmer(data = PM_G1, BetaValue ~ meth.diff.parentals : Treatment + (1|CpGSite) + (1|brotherPairID))

## check normality of residuals assumption
qqnorm(resid(parmod))
qqline(resid(parmod))

pred <- ggpredict(parmod, terms = c("meth.diff.parentals", "Treatment"))
plot(pred, add.data = T)+
  scale_color_manual(values = c("black", "red"))+
  scale_y_continuous(name = "Beta values")+
  scale_x_continuous(name = "Methylation difference between infected and control parents in percentage")+
  ggtitle("Predicted methylation ratio (Beta) values in parents\n as a function of differential methylation between exposed and control groups")+
  theme_bw()

##############
## Linear model: does the beta value of offspring at DMS depends on treatment Parent x Offspring?
##############
modFull <- lmer(BetaValue ~ (G1_trt * G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID),data = PM_G2, REML = F) # REML =F for model comparison
mod_noG1trt <- lmer(BetaValue ~ G2_trt:hypohyper + (1|CpGSite)+ (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
mod_noG2trt <-lmer(BetaValue ~ G1_trt:hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
mod_noInteractions <- lmer(BetaValue ~ (G1_trt + G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)
mod_noHypoHyper <- lmer(BetaValue ~ (G1_trt * G2_trt) + (1|CpGSite) + (1|Sex) + (1|brotherPairID), data = PM_G2, REML = F)

## check normality of residuals assumption
qqnorm(resid(modFull))
qqline(resid(modFull))

## Likelihood ratio tests for all variables:
lrtest(modFull, mod_noG1trt) # G1 trt is VERY VERY significant (LRT: χ² (4) = 1163.6, p < 0.001)
lrtest(modFull, mod_noG2trt) # G2 trt is VERY VERY significant (LRT: χ² (4) = 30.02, p < 0.001) NB that changed when brotherpair is used instead of family!
lrtest(modFull, mod_noInteractions) # interactions are significant (LRT: χ² (2) = 9.21, p < 0.01)
lrtest(modFull, mod_noHypoHyper) # hypo/hyper VERY VERY significant (LRT: χ² (4) = 1140, p < 0.001)

## Post-hoc tests between treatments
modFull <- lmer(BetaValue ~ (G1_trt * G2_trt):hypohyper + (1|CpGSite) + (1|Sex) + (1|brotherPairID),data = PM_G2)
modFull_emmeans <- emmeans(modFull, list(pairwise ~ (G1_trt:G2_trt):hypohyper), adjust = "tukey")
modFull_emmeans

P1 <- plot(modFull_emmeans, by = "hypohyper", comparisons = TRUE) +
  # coord_flip()+
  theme_bw() +
  ggtitle("Estimated marginal means of methylation ratio (beta)\n of offspring at parental DMS")+
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_continuous("Beta value (methylation ratio)", limits = c(47,69.5))

## NB: Comparison arrows: https://cran.r-project.org/web/packages/emmeans/vignettes/xplanations.html
## two estimated marginal means (EMMs) differ significantly if, and only if, their respective comparison arrows do not overlap
## These comparison arrows are decidedly not the same as confidence intervals.
## Confidence intervals for EMMs are based on the statistical properties of the individual EMMs, whereas comparison arrows
## are based on the statistical properties of differences of EMMs.

## Add the PARENTAL DMS value
## Same test on ALL, G1 and G2 fish
modFullG1 <- lmer(BetaValue ~ G1_trt:hypohyper + (1|CpGSite) + (1|brotherPairID), data = PM_G1)

modFullG1_emmeans <- emmeans(modFullG1, list(pairwise ~ G1_trt:hypohyper), adjust = "tukey")
modFullG1_emmeans

P2 <- plot(modFullG1_emmeans, by = "hypohyper", comparisons = TRUE) +
  theme_bw() +
  ggtitle("Estimated marginal means of methylation ratio (beta)\n of parents at DMS")+
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_continuous("Beta value (methylation ratio)", limits = c(47,69.5))

ggarrange(P2, P1, labels = c("A", "B"), ncol = 1, nrow = 2)

#####################################################################################
## Are the nbr of residuals methylation AT PARENTAL DMS different in the 4 G2 trt? ##
## (for hypo vs hypermeth)? ##
##############################

length(unique(PM_G1$CpGSite))# 3648 positions
PM_G1 %>% dplyr::count(ID)## NB: not all covered in all samples

length(unique(PM_G2$CpGSite[PM_G2$hypohyper %in% "hypo"]))# 1176 positions hypomethylated upon parental inf
length(unique(PM_G2$CpGSite[PM_G2$hypohyper %in% "hyper"]))# 2472 positions hypermethylated upon parental inf

myfun <- function(X){
  ## Calculate nbr of CpG hypo/hypermethylated per individual, and nbr of covered CpG:
  X <- X %>% group_by(ID, Treatment, brotherPairID, clutch.ID, Sex) %>%
    dplyr::summarise(ncov = n(),
                     hypoMeth = sum(BetaValue < 0.3),
                     hyperMeth = sum(BetaValue > 0.7)) %>% data.frame()
  # Calculate residuals of nbr of methCpG by nbr of covered CpG
  X$res_Nbr_methCpG_Nbr_coveredCpG_HYPO = residuals(lm(X$hypoMeth ~ X$ncov))
  X$res_Nbr_methCpG_Nbr_coveredCpG_HYPER = residuals(lm(X$hyperMeth ~ X$ncov))

  ## Statistical model: is it different in each offspring trt group?
  mod1 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
               data = X, REML = F)
  mod0 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ 1 + (1|brotherPairID/clutch.ID) + (1|Sex),
               data = X, REML = F)
  print(lrtest(mod1, mod0))

  ## Post-hoc test:
  modhypo <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPO ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
                  data = X)
  ## pairwise posthoc test
  print(emmeans(modhypo, list(pairwise ~ Treatment), adjust = "tukey"))

  mod3 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
               data = X, REML = F)
  mod4 <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ 1 + (1|brotherPairID/clutch.ID) + (1|Sex),
               data = X, REML = F)
  print(lrtest(mod3, mod4))

  ## Post-hoc test:
  modhyper <- lmer(res_Nbr_methCpG_Nbr_coveredCpG_HYPER ~ Treatment + (1|brotherPairID/clutch.ID) + (1|Sex),
                   data = X)
  ## pairwise posthoc test
  print(emmeans(modhyper, list(pairwise ~ Treatment), adjust = "tukey"))

  ## Plot
  phypo <- plot(ggpredict(modhypo, terms = c("Treatment")), add.data = TRUE)+
    scale_y_continuous("Residuals of number of hypomethylated methylated \ncytosines on number of cytosines covered") +
    ggtitle("Predicted residuals nbr of hypomethylated CpG")+
    theme_bw()

  phyper <- plot(ggpredict(modhyper, terms = c("Treatment")), add.data = TRUE)+
    scale_y_continuous("Residuals of number of hypermethylated methylated \n cytosines on number of cytosines covered") +
    ggtitle("Predicted residuals nbr of hypermethylated CpG")+
    theme_bw()
  return(list(phypo, phyper))
}

listplots <- myfun(X = PM_G2[PM_G2$hypohyper %in% "hypo",])
## NOT significant
annotate_figure(ggarrange(listplots[[1]], listplots[[2]],ncol = 2, nrow = 1),
                top = text_grob("Parental DMS hypomethylated upon infection, in offspring"))

listplots <- myfun(X = PM_G2[PM_G2$hypohyper %in% "hyper",])
## Treatment SIGNIFICANT in both excess hypo/hyper methylation **

# $`pairwise differences of Treatment`
## HYPO
# 1                       estimate   SE   df t.ratio p.value
# NE_control - E_control     23.71 7.28 10.3   3.257  0.0353
# NE_control - E_exposed     26.88 7.26 10.3   3.701  0.0172

## HYPER
# 1                       estimate   SE   df t.ratio p.value
# NE_control - E_control    -24.06 7.36 10.3  -3.269  0.0348
# NE_control - E_exposed    -27.07 7.34 10.3  -3.687  0.0177

annotate_figure(ggarrange(listplots[[1]], listplots[[2]],ncol = 2, nrow = 1),
                top = text_grob("Parental DMS hypermethylated upon infection, in offspring"))

## --> The beta values in parentalDMS in offspring follow the parental pattern hypo/hyper methylated upon infection

# plot beta and parent 

## How to find WHICH POSITIONS ARE TRANSMISSIBLE STATES???

## BUG
# PMbyclutch <- merge(data.frame(clutch.ID = PM_G2$clutch.ID, betaG2 = PM_G2$BetaValue, G2_trt = PM_G2$G2_trt),
#                     data.frame(clutch.ID = PM_G1$clutch.ID, betaG1 = PM_G1$BetaValue, G1_trt = PM_G1$G1_trt))

##############
## Plot mean Beta values of offsprings at parental DMS, per trt, along the chromosomes
##############
meanBeta_G2_simple <- PM_G2 %>% group_by(Chr, Pos, Treatment) %>%
  dplyr::summarize(Mean = mean(BetaValue, na.rm=TRUE))
names(meanBeta_G2_simple) <- c("Chr", "Pos", "Treatment_G2", "MeanBetaG2")

genome <- GYgynogff %>%
  mutate(chrom_nr=chrom %>% deroman(), chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
  arrange(chrom_order) %>%
  mutate(gstart=lag(length,default=0) %>% cumsum(), gend=gstart+length, 
         type=LETTERS[2-(chrom_order%%2)],   gmid=(gstart+gend)/2)

mydata = tibble(trt=meanBeta_G2_simple$Treatment_G2,
                chrom=meanBeta_G2_simple$Chr, pos=meanBeta_G2_simple$Pos, beta=meanBeta_G2_simple$MeanBetaG2)
mydata$pos <- as.numeric(mydata$pos)
mydata$pos <- as.numeric(mydata$pos)

table(mydata$chrom)## check that chrXIX and chrUN are well removed!!

# join DMS and genomic position
data = dplyr::left_join(mydata, genome) %>% dplyr::mutate(gpos=pos+gstart)

# plot:
ggplot()+
  geom_rect(data=genome,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
  geom_point(data=data, aes(x=gpos,y=beta, shape=trt, col=trt),fill="white", size = 1)+
  scale_color_manual(values = colOffs)+
  scale_shape_manual(values=c(21,21,21,21))+
  scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
  scale_x_continuous(breaks=genome$gmid,labels=genome$chrom %>% str_remove(.,"Gy_chr"),
                     position = "top",expand = c(0,0))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line=element_blank(),
        axis.title = element_blank(),
        strip.placement = "outside")+
  facet_grid(trt~.)+
  ggtitle("Mean methylation proportions at the 3648 parental DMS for each offspring group")

###########################################################################################
## Not conclusive: test hyp that beta value is different (higher?) in the center of the chromosome, 
## where there are less recombinations. Test for each chromosome
meanBeta_G2_extended <- PM_G2 %>% group_by(Chr, Pos, Treatment, Sex, brotherPairID) %>%
  dplyr::summarize(Mean = mean(BetaValue, na.rm=TRUE))
names(meanBeta_G2_extended) <- c("Chr", "Pos", "Treatment_G2","Sex", "brotherPairID", "MeanBetaG2")

mydata = tibble(trt=meanBeta_G2_extended$Treatment_G2, sex = meanBeta_G2_extended$Sex, brotherPairID = meanBeta_G2_extended$brotherPairID,
                chrom=meanBeta_G2_extended$Chr, pos=meanBeta_G2_extended$Pos, beta=meanBeta_G2_extended$MeanBetaG2)
mydata$pos <- as.numeric(mydata$pos)

# join DMS and genomic position
data = dplyr::left_join(mydata, genome) %>% dplyr::mutate(gpos=pos+gstart)

## Add distance to center
data$dist2mid <- abs(data$gmid - data$gpos)

mod <- lmer(beta ~ dist2mid:chrom + (1|trt) + (1|sex) + (1|brotherPairID), data = data, REML = F)
mod0 <- lmer(beta ~ dist2mid + (1|trt) + (1|sex) + (1|brotherPairID), data = data, REML = F)
mod00 <- lmer(beta ~ chrom + (1|trt) + (1|sex) + (1|brotherPairID), data = data, REML = F)

lrtest(mod, mod0) # chromosome matters
lrtest(mod0, mod00) # distance to middle matters

## check normality of residuals assumption
qqnorm(resid(mod))
qqline(resid(mod)) # quite skewed

pred <- ggpredict(mod, terms = c("dist2mid","chrom"))
plot(pred) +
  scale_color_manual(values = 1:20)

#########################
## Clustering analysis ##
#########################
## Run Adonis to this if clustering is done by treatment
x = percMethylation(uniteCov14_G2_woSexAndUnknowChrOVERLAP)

## Transpose
x=t(x)

## Creates a distance matrix. Method: Bray-Curtis, package vegan
data.dist = as.matrix(vegdist(x, "bray", upper = FALSE, na.rm = T))

## Check that the order is the same than with the metadata
table(fullMetadata_OFFS$SampleID == rownames(data.dist))

# We use a PERMANOVA to test the hypothesis that paternal treatment, 
# offspring treatment, sex and their interactions significantly influencing global methylation
perm <- how(nperm = 1000) # 1000 permutations
setBlocks(perm) <- with(fullMetadata_OFFS, brotherPairID) # define the permutation structure considering brotherPairID

## Full model
adonis2(data.dist ~ PAT * outcome * Sex, data = fullMetadata_OFFS, permutations = perm)
# adonis2(formula = data.dist ~ PAT * outcome * Sex, data = fullMetadata_OFFS, permutations = perm)
#                   Df SumOfSqs      R2      F   Pr(>F)    
# PAT               1 0.004138 0.01330 1.4748 0.000999 ***
# outcome           1 0.002950 0.00948 1.0515 0.075924 .  
# Sex               1 0.003755 0.01207 1.3383 0.004995 ** 
# PAT:outcome       1 0.002869 0.00922 1.0227 0.050949 .  

##### NMDS
# find the best number of dimensions (goeveg lib)
## Clarke 1993 suggests the following guidelines for acceptable stress values: <0.05 = excellent, <0.10
# = good, <0.20 = usable, >0.20 = not acceptable. The plot shows the border of the 0.20 stress value
# limit. Solutions with higher stress values should be interpreted with caution and those with stress
# above 0.30 are highly suspect
dimcheckMDS(
  data.dist,
  distance = "bray",
  k = 7,
  trymax = 100,
  autotransform = TRUE
)
abline(h = 0.1, col = "darkgreen")

#Create NMDS based on bray-curtis distances - metaMDS finds the
# most stable NMDS solution by randomly starting from different points in your data
set.seed(1234)

# generate CpG sums these will be NA is any are missing
csum <- colSums(x)
# check if any are missing
any(is.na(csum))
# TRUE
# yes, some missing, so which ones?
length(which(!is.na(csum)))#78384 complete positions

x = x[,which(!is.na(csum))]

NMDS <- metaMDS(comm = x, distance = "bray", maxit=1000, k = 3)

#check to see stress of NMDS
mystressplot <- stressplot(NMDS) 

#extract plotting coordinates
MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]
## OR #extract NMDS scores (x and y coordinates)
## data.scores = as.data.frame(scores(NMDS))

#create new data table (important for later hulls finding)
# with plotting coordinates and variables to test (dim 1,2,3)
NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                 trtG1G2 = fullMetadata_OFFS$trtG1G2)

# generating convex hulls splitted by myvar in my metadata:
hulls <- NMDS_dt[, .SD[chull(MDS1, MDS2)], by = trtG1G2]

## to insert image (expe design)
img <- readPNG("../../data/designExpeSimple.png")
g <- rasterGrob(img, interpolate=TRUE)

myNMDSplot <- ggplot(NMDS_dt, aes(x=MDS1, y=MDS2)) +
  geom_polygon(data = hulls, aes(fill=trtG1G2), alpha=0.3) +
  # scale_color_manual(values = colOffs)+
  scale_fill_manual(values = colOffs)+
  geom_point(aes(fill=trtG1G2, shape=trtG1G2), col = "black", size = 3, alpha =.5) +
  scale_shape_manual(values = c(21,22,23,24)) +
  # geom_label(aes(label=rownames(NMDS$points), col = fullMetadata_OFFS$brotherPairID))+
  theme_bw() +
  # theme(legend.position = "none") +
  annotation_custom(g, xmin=-0.25, xmax=-0.1, ymin=0.05, ymax=0.15) +
  ggtitle("NMDS analysis of the 78384 CpG sites from parental control/infected DMS,\npresent in all offspring")
myNMDSplot

################################################
################ Venn diagrams #################
################################################
## Just considering intersecting CpGs:
myVennFUN <- function(A, B, C, catnames, myCols = c("grey","green","red")){
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")#to rm log files
  Venn <- venn.diagram(
    x = list(A, B, C), category.names = catnames, filename = NULL,
    margin = 0, lwd = 2, lty = 'blank', fill = myCols,
    cex = .4, fontface = "bold",fontfamily = "sans", print.mode=c("raw","percent"),
    cat.cex = 0.4, cat.fontface = "bold", cat.default.pos = "outer",
    cat.col = myCols, cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.055),
    cat.fontfamily = "sans", rotation = 1
  )
  return(Venn)
}

## output Venn diagrams
allVenn <- myVennFUN(A = DMS_info_G1$DMS, B = DMS_info_G2_G1c_final$DMS, C = DMS_info_G2_G1i_final$DMS,
                     catnames = c("DMS G1" , "DMS G2-c", "DMS G2-i"))
hypoVenn <- myVennFUN(A = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hypo"],
                      B = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hypo"],
                      C = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hypo"],
                      catnames = c("DMS G1\nhypo" , "DMS G2-c\nhypo", "DMS G2-i\nhypo"))
hyperVenn <- myVennFUN(A = DMS_info_G1$DMS[DMS_info_G1$direction %in% "hyper"],
                       B = DMS_info_G2_G1c_final$DMS[DMS_info_G2_G1c_final$direction %in% "hyper"],
                       C = DMS_info_G2_G1i_final$DMS[DMS_info_G2_G1i_final$direction %in% "hyper"],
                       catnames = c("DMS G1\nhyper" , "DMS G2-c\nhyper", "DMS G2-i\nhyper"))

# Output the diagrams
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_hyper_hypo.png", width = 5.5, height = 4.5, units = 'in', res = 300)
pushViewport(plotViewport(layout=grid.layout(2, 3)))
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(allVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(hypoVenn)
popViewport()
pushViewport(plotViewport(layout.pos.col=3, layout.pos.row=2))
grid.draw(hyperVenn)
dev.off()

rm(allVenn, hypoVenn, hyperVenn)

######################
## Features Annotation (use package genomation v1.24.0)
## NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
## The default option is to take -1000,+1000bp around the TSS and you can change that. 
## -> following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream

#########################
## load genome annotation
gene.obj=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.bed12",
                                remove.unusual = FALSE, up.flank = 1500, down.flank = 500)

gene.obj.2=readTranscriptFeatures("../../gitignore/bigdata/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_afterAGATcuration.bed12",
                                  remove.unusual = FALSE, up.flank = 1500, down.flank = 500)


par(mfrow=c(1,3))
par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
## Parents comparison:
diffAnn_PAR = annotateWithGeneParts(as(DMS15pc_G1_half,"GRanges"),gene.obj)
diffAnn_PAR
plotTargetAnnotation(diffAnn_PAR,precedence=TRUE, main="DMS G1", 
                     cex.legend = 1, border="white")

## Offspring from control parents comparison:
diffAnn_G2_controlG1 = annotateWithGeneParts(as(DMS15pc_G2_controlG1_half,"GRanges"),gene.obj)
diffAnn_G2_controlG1
plotTargetAnnotation(diffAnn_G2_controlG1,precedence=TRUE, main="DMS G2-G1c", 
                     cex.legend = 1, border="white")
## Offspring from infected parents comparison:
diffAnn_G2_infectedG1 = annotateWithGeneParts(as(DMS15pc_G2_infectedG1_half,"GRanges"),gene.obj)
diffAnn_G2_infectedG1
plotTargetAnnotation(diffAnn_G2_infectedG1,precedence=TRUE, main="DMS G2-G1i", 
                     cex.legend = 1, border="white")
par(mfrow=c(1,1)) 

##########################
## Separate hyper and hypo
runHyperHypoAnnot <- function(){
  par(mfrow=c(2,3))
  par(mar = c(.1,0.1,5,0.1)) # Set the margin on all sides to 2
  ####### HYPO
  ## Parents comparison:
  A = annotateWithGeneParts(
    as(DMS15pc_G1_half[DMS_info_G1$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(A,precedence=TRUE, main="DMS G1\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  B = annotateWithGeneParts(
    as(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(B,precedence=TRUE, main="DMS G2-G1c\nhypo", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  C = annotateWithGeneParts(
    as(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hypo",],"GRanges"),gene.obj)
  plotTargetAnnotation(C,precedence=TRUE, main="DMS G2-G1i\nhypo", 
                       cex.legend = .4, border="white")
  ####### HYPER
  ## Parents comparison:
  D = annotateWithGeneParts(
    as(DMS15pc_G1_half[DMS_info_G1$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(D,precedence=TRUE, main="DMS G1\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from control parents comparison:
  E = annotateWithGeneParts(
    as(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(E,precedence=TRUE, main="DMS G2-G1c\nhyper", 
                       cex.legend = .4, border="white")
  ## Offspring from infected parents comparison:
  f = annotateWithGeneParts(
    as(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hyper",],"GRanges"),gene.obj)
  plotTargetAnnotation(f,precedence=TRUE, main="DMS G2-G1i\nhyper", 
                       cex.legend = .4, border="white")
  par(mfrow=c(1,1))
  return(list(G1hypo=A, G2G1chypo=B, G2G1ihypo=C, G1hyper=D, G2G1chyper=E, G2G1ihyper=f))
}

myannot=runHyperHypoAnnot()

############################################################
## Venn diagram of overlapping features by their annotation:
table(rowSums(as.data.frame(myannot$G1hypo@members))) # NB: some positions are labelled with several features!
## as in MBE 2021: "giving precedence to the following order promoters, exons,
## introns, and intergenic regions when features overlapped"

myAnnotateDMS <- function(DMS, annot){
  ## sanity check
  if (nrow(DMS) != nrow(annot)){"STOP error in arguments"}
  DMS$pos <- paste(DMS$chr, DMS$start, DMS$end)
  ## NB as in MBE 2021: "giving precedence to the following order promoters, exons,
  ## introns, and intergenic regions when features overlapped"
  DMS$feature <- NA
  ## 1. promoters
  DMS$feature[which(annot$prom == 1)] = "promoter"
  ## 2. exons
  DMS$feature[which(annot$exon == 1 & annot$prom ==0)] = "exon"
  ## 3. intron
  DMS$feature[which(annot$intro == 1 & annot$exon == 0 & annot$prom ==0)] = "intron"
  ## 4. intergenic regions
  DMS$feature[which(annot$intro == 0 & annot$exon == 0 & annot$prom ==0)] = "intergenic"
  return(DMS)
}

DMS15pc_G1_half = myAnnotateDMS(DMS15pc_G1_half, as.data.frame(diffAnn_PAR@members))
DMS15pc_G1_half_HYPO = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hypo",],
                                     as.data.frame(myannot$G1hypo@members))
DMS15pc_G1_half_HYPER = myAnnotateDMS(DMS15pc_G1_half[DMS_info_G1$direction %in% "hyper",],
                                      as.data.frame(myannot$G1hyper@members))

DMS15pc_G2_controlG1_half = myAnnotateDMS(DMS15pc_G2_controlG1_half, as.data.frame(diffAnn_G2_controlG1@members))
DMS15pc_G2_controlG1_half_HYPO = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hypo",],
                                               as.data.frame(myannot$G2G1chypo@members))
DMS15pc_G2_controlG1_half_HYPER = myAnnotateDMS(DMS15pc_G2_controlG1_half[DMS_info_G2_G1c_final$direction %in% "hyper",],
                                                as.data.frame(myannot$G2G1chyper@members))

DMS15pc_G2_infectedG1_half = myAnnotateDMS(DMS15pc_G2_infectedG1_half, as.data.frame(diffAnn_G2_infectedG1@members))
DMS15pc_G2_infectedG1_half_HYPO = myAnnotateDMS(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hypo",],
                                                as.data.frame(myannot$G2G1ihypo@members))
DMS15pc_G2_infectedG1_half_HYPER = myAnnotateDMS(DMS15pc_G2_infectedG1_half[DMS_info_G2_G1i_final$direction %in% "hyper",],
                                                 as.data.frame(myannot$G2G1ihyper@members))

## Make Venn diagram for each feature
getFeatureDFHYPO <- function(myfeat){
  a = DMS15pc_G1_half_HYPO$pos[DMS15pc_G1_half_HYPO$feature %in% myfeat]
  b = DMS15pc_G2_controlG1_half_HYPO$pos[DMS15pc_G2_controlG1_half_HYPO$feature %in% myfeat]
  c = DMS15pc_G2_infectedG1_half_HYPO$pos[DMS15pc_G2_infectedG1_half_HYPO$feature %in% myfeat]
  return(list(a=a,b=b,c=c))
}

getFeatureDFHYPER <- function(myfeat){
  a = DMS15pc_G1_half_HYPER$pos[DMS15pc_G1_half_HYPER$feature %in% myfeat]
  b = DMS15pc_G2_controlG1_half_HYPER$pos[DMS15pc_G2_controlG1_half_HYPER$feature %in% myfeat]
  c = DMS15pc_G2_infectedG1_half_HYPER$pos[DMS15pc_G2_infectedG1_half_HYPER$feature %in% myfeat]
  return(list(a=a,b=b,c=c))
}

getVenn <- function(feat, direction){
  if (direction == "hypo"){
    myVennFUN(getFeatureDFHYPO(feat)[["a"]],
              getFeatureDFHYPO(feat)[["b"]],
              getFeatureDFHYPO(feat)[["c"]],
              catnames = c(paste0("DMS G1\nhypo\n", feat), 
                           paste0("DMS G2-c\nhypo\n", feat),
                           paste0("DMS G2-i\nhypo\n", feat)))
  }else if (direction == "hyper"){
    myVennFUN(getFeatureDFHYPER(feat)[["a"]],
              getFeatureDFHYPER(feat)[["b"]],
              getFeatureDFHYPER(feat)[["c"]],
              catnames = c(paste0("DMS G1\nhyper\n", feat), 
                           paste0("DMS G2-c\nhyper\n", feat),
                           paste0("DMS G2-i\nhyper\n", feat)))
  }
}

## output Venn diagrams
png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_byfeatures_HYPO.png", width = 5, height = 5.5, units = 'in', res = 300)
#height = 500, width = 500, compression = "lzw")
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(getVenn("promoter", "hypo"))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(getVenn("exon", "hypo"))
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intron", "hypo")))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intergenic", "hypo"))) 
dev.off()

png(file="Rfigures/VennDMSinhalffish_intersectingCpGs_byfeatures_HYPER.png", width = 5, height = 5.5, units = 'in', res = 300)
#height = 500, width = 500, compression = "lzw")
pushViewport(plotViewport(layout=grid.layout(2, 2)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(getVenn("promoter", "hyper"))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=1))
grid.draw(getVenn("exon", "hyper"))
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intron", "hyper")))
popViewport()
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=2))
grid.draw(grid.draw(getVenn("intergenic", "hyper"))) 
dev.off()

###########################
## Manhattan plot of DMS ##
###########################
## Parents trt-ctrl
# load annotation
annot_PAR <- as.data.frame(diffAnn_PAR@members)
makeManhattanPlots(DMSfile = DMS15pc_G1_half, annotFile = annot_PAR, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

## G2-G1c trt-ctrl
# load annotation
annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_controlG1_half, annotFile = annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

## G2-G1i trt-ctrl
# load annotation
annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)
makeManhattanPlots(DMSfile = DMS15pc_G2_infectedG1_half, annotFile = annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

## Outliers in Manhattan plot: 15% diff + 2SD
outliers_G1_final <- which(abs(DMS15pc_G1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G1_half$meth.diff)))
outliers_annot_G1 <- as.data.frame(diffAnn_PAR@members)[outliers_G1_final,]
makeManhattanPlots(DMSfile = DMS15pc_G1_half[outliers_G1_final, ],
                   annotFile = outliers_annot_G1, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G1 DMS")

outliers_G2_G1c_final <- which(abs(DMS15pc_G2_controlG1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G2_controlG1_half$meth.diff)))
outliers_annot_G2_G1c <- as.data.frame(diffAnn_G2_controlG1@members)[outliers_G2_G1c_final,]
makeManhattanPlots(DMSfile = DMS15pc_G2_controlG1_half[outliers_G2_G1c_final, ],
                   annotFile = outliers_annot_G2_G1c, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1c DMS")

outliers_G2_G1i_final <- which(abs(DMS15pc_G2_infectedG1_half$meth.diff) > 15 + 2*sd(abs(DMS15pc_G2_infectedG1_half$meth.diff)))
outliers_annot_G2_G1i <- as.data.frame(diffAnn_G2_infectedG1@members)[outliers_G2_G1i_final,]
makeManhattanPlots(DMSfile = DMS15pc_G2_infectedG1_half[outliers_G2_G1i_final, ],
                   annotFile = outliers_annot_G2_G1i, GYgynogff = GYgynogff, 
                   mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of G2-G1i DMS")

###############################
## Gene association with DMS ##
###############################

## Keep CpG on genes:

## Heckwolf 2020: To be associated to a gene, the pop-DMS had to be either inside the gene or, 
## if intergenic, not further than 10 kb away from the TSS. 
rows2rm <- which(diffAnn_PAR@dist.to.TSS$dist.to.feature>10000 & 
                   rowSums(diffAnn_PAR@members) %in% 0)

DMS15pc_G1_half_GENES <- DMS15pc_G1_half[-rows2rm,]

# 4393 to keep, 681 to rm
############ Follow up
gene.obj$promoters

DMS15pc_G1_half_GENES





setwd("/home/user/R/x86_64-pc-linux-gnu-library/4.0/genomation/extdata")
gff.file=system.file("extdata/Bubalus_bubalis.gtf", package = "genomation")
gff = gffToGRanges(gff.file)
head(gff)
Anno <- GRangesList(gff)
setwd("~/Desktop/Final_sort")
diffAnnhyper3=annotateWithGeneParts(as(myDiff25p.hyper,"Granges"),Anno)


?readTranscriptFeatures
head(getAssociationWithTSS(diffAnn_PAR))


### BONUS
# ##############
# ## Beta values of offsprings at parental DMS, per trt, along the percent DM
# ##############
# dfHypo <- PM_G2[PM_G2$hypohyper %in% "hypo",]
# dfHyper <- PM_G2[PM_G2$hypohyper %in% "hyper",]
# 
# ### PLOTS
# ## to insert image (expe design)
# img <- readPNG("../../data/designExpeSimple.png")
# g <- rasterGrob(img, interpolate=TRUE)
# 
# ## Hypo
# mod <- lmer(BetaValue ~ meth.diff.parentals:Treatment + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#             data = dfHypo, REML = F)
# mod0 <- lmer(BetaValue ~ meth.diff.parentals + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#              data = dfHypo, REML = F)
# lrtest(mod, mod0) # trt significant chisq=549.14  p < 2.2e-16 ***
# 
# mod <- lmer(BetaValue ~ meth.diff.parentals:Treatment + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#             data = dfHypo, REML = T)
# predPosBeta <- ggpredict(mod, terms = c("meth.diff.parentals", "Treatment"))
# plotHypo <- plot(predPosBeta)+
#   scale_color_manual(values = colOffs) +
#   theme_bw()+
#   theme(legend.position = "none")+
#   scale_x_continuous(limits = c(-60,-15))+
#   scale_y_continuous(limits = c(35,65))+
#   ggtitle("Predicted methylation percentage in offspring at parental hypomethylated DMS")
# 
# ## Hyper
# mod <- lmer(BetaValue ~ meth.diff.parentals:Treatment + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#             data = dfHyper, REML = F)
# mod0 <- lmer(BetaValue ~ meth.diff.parentals + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#              data = dfHyper, REML = F)
# lrtest(mod, mod0) # trt significant chisq=662.12, p < 2.2e-16 ***
# 
# mod <- lmer(BetaValue ~ meth.diff.parentals:Treatment + (1|CpGSite) + (1|Sex) + (1|brotherPairID), 
#             data = dfHyper, REML = T)
# predPosBeta <- ggpredict(mod, terms = c("meth.diff.parentals", "Treatment"))
# plotHyper <- plot(predPosBeta)+
#   scale_color_manual(values = colOffs) +
#   theme_bw()+
#   theme(legend.position = "none") +
#   annotation_custom(g, xmin=35, xmax=60, ymin=55, ymax=65)+
#   scale_x_continuous(limits = c(15,60))+
#   scale_y_continuous(limits = c(35,65))+
#   ggtitle("Predicted methylation percentage in offspring at parental hypermethylated DMS")
# 
# grid.arrange(plotHypo, plotHyper, ncol=2)