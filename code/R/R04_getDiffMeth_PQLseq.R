## Data preparation for PQLseq
## A. Balard
## 2nd of August 2024

## load previous scripts
source("R03_prepObjectMethylkit_runInCLUSTER.R")

######################################
## load objects created in previous scripts
## Parents file filtered and normalised with methylkit
load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovHALF_G1_woSexAndUnknowChr_OVERLAPwG2_20dec2022.RData")

## Sanity check
uniteCovHALF_G1_woSexAndUnknowChrOVERLAP %>% nrow == 1002565

## Offspring file filtered and normalised with methylkit
load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovHALF_G2_woSexAndUnknowChr_OVERLAPwG1_20dec2022.RData")

## Sanity check
uniteCovHALF_G2_woSexAndUnknowChrOVERLAP %>% nrow == 1002565

##################################
## PQLseq differential analysis ##
##################################

## PQLseq requires four input files containing (1) methylated read counts, (2) total read counts, (3) relatedness matrix, and (4) predictor variable of interest

######################
## relatedness matrix: 1 against oneself, 0.5 between sibling, 0.25 uncle/nephew-niece, 0.125 cousins
fullMetadata = read.csv("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/data/fullMetadata127_Alice.csv")

# order by father family, brother pair then clutch
met = fullMetadata
met = met[order(met$FamilyOfFather, met$brotherPairID, met$clutch.ID),]

relatMat <- matrix(nrow=nrow(met), ncol=nrow(met), data = 0)
rownames(relatMat)<-met$SampleID
colnames(relatMat)<-met$SampleID

## coef of relatedness with oneself = 1
diag(relatMat) <- 1

## Same clutch: coef of relatedness with your sibling/father/offspring = 0.5 
for (i in 1:nrow(met)) {
  for (j in 1:nrow(met)) {
    if (relatMat[i, j]==0 && met$clutch.ID[i] == met$clutch.ID[j]) {
      relatMat[i, j] = 0.5
    }}
}
## Same BP, diff clutch, gener P: coef of relatedness with your brother in G1 = 0.5
for (i in 1:nrow(met)) {
  for (j in 1:nrow(met)) {
    if (relatMat[i, j]==0 &&  met$brotherPairID[i] == met$brotherPairID[j] &&
        met$Generat[i] == met$Generat[j] && met$Generat[i] == "P") {
      relatMat[i, j] = 0.5
    }}
}
## Same BP, diff clutch, gener O: coef of relatedness with your cousins in G2 = 0.125
for (i in 1:nrow(met)) {
  for (j in 1:nrow(met)) {
    if (relatMat[i, j] == 0 && met$brotherPairID[i] == met$brotherPairID[j] &&
        met$Generat[i] == met$Generat[j] && met$Generat[i] == "O"){
      relatMat[i, j] = 0.125
    }}
}
## Same BP, diff clutch, diff gener: coef of relatedness with your uncle/niece/nephew = 0.25
for (i in 1:nrow(met)) {
  for (j in 1:nrow(met)) {
    if (relatMat[i, j] == 0 && met$brotherPairID[i] == met$brotherPairID[j]&&
        met$Generat[i] != met$Generat[j]) {
      relatMat[i, j] = 0.25
    }}
}


## Write out
# write.csv( relatMat

## Sanity check (to do in IDE)


## load ggcorr function instead of full GGally
source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R")
ggcorr(data = NULL, cor_matrix = relatMat, size = 2, color = "grey50")
## Ok, all the relationships within the 8 families are well represented





makefiles.1.2 <- function(uniteCov){
    ## C represents a methylated cytosine that remained unchanged during bisulfite treatment
    dfC=getData(uniteCov)[ ,grep("numCs", names(uniteCov))]
    colnames(dfC)=uniteCov@sample.ids

    dftot=getData(uniteCov)[ ,grep("coverage", names(uniteCov))]
    colnames(dftot)=uniteCov@sample.ids

    counts=cbind(data.frame(siteID = paste(uniteCov$chr, uniteCov$start, sep ="_")), dfC)
    totalcounts=cbind(data.frame(siteID = paste(uniteCov$chr, uniteCov$start, sep ="_")), dftot)

    return(list(counts=counts, totalcounts=totalcounts))
}

G1 <- makefiles.1.2(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP)
counts_G1 = G1$counts
totalcounts_G1 = G1$totalcounts

G2 <- makefiles.1.2(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP)
counts_G2 = G2$counts
totalcounts_G2 = G2$totalcounts

### 5 comparisons to do:
# 0. PARENTS trt-ctrl
# 1. CC-TC = CONTROL fish (parent CvsT)
# 2. CT-TT = TREATMENT fish (parent CvsT)
# 3. CC-CT = fish from CONTROL parents (G2 CvsT)
# 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

## 1&2: paternal effect 
## 3&4: G2 direct treatment effect

################################################################
## Correspondance trtG1G2 and numerical values used by MethylKit
# table(fullMetadata$trtG1G2_NUM, fullMetadata$trtG1G2)
## Control C 1
## Exposed T 4
## NE_control CC 5
## NE_exposed CT 6
## E_control TC 2
## E_exposed TT 3

#####################################################################
## Calculate DMS/DMR accounting for covariates: brotherPairID and sex

## (4) phenotype file

## parents: 
uniteCovHALF_G1_woSexAndUnknowChrOVERLAP@treatment

## offspring:
uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@treatment

## select the correct counts, totalcounts, and phenotype for a given comparison
makesubcountsG2 <- function(groupa, groupb){

    subpheno = uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@treatment[uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@treatment %in% c(groupa,groupb)]

    subcounts = counts_G2[names(counts_G2) %in% c("siteID", uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@sample.ids[uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@treatment %in% c(groupa,groupb)])]
    
    subtotalcounts = totalcounts_G2[names(totalcounts_G2) %in% c("siteID", uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@sample.ids[uniteCovHALF_G2_woSexAndUnknowChrOVERLAP@treatment %in% c(groupa,groupb)])]
    return(list(comp=paste("Comparison between group ", groupa, " and ", groupb),
                subpheno=subpheno, subcounts=subcounts, subtotalcounts=subtotalcounts))
}

## 1. CC-TC = CONTROL fish (parent CvsT)
G2_CC.TC = makesubcountsG2(2,5)

G2_CC.TC$subpheno
G2_CC.TC$subcounts
G2_CC.TC$subtotalcounts

## 2. CT-TT = TREATMENT fish (parent CvsT)

## 3. CC-CT = fish from CONTROL parents (G2 CvsT)

## 4. TC-TT = fish from TREATMENT parents (G2 CvsT)

