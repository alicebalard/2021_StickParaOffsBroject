## Data preparation for PQLseq
## A. Balard
## 2nd of August 2024

## libraries
library(PQLseq)
library(dplyr)

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

uniteCovHALF_G2_woSexAndUnknowChrOVERLAP %>% head

##################################
## PQLseq differential analysis ##
##################################

## PQLseq requires four input files containing (1) methylated read counts, (2) total read counts, (3) relatedness matrix, and (4) predictor variable of interest

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
uniteCovHALF_G1_woSexAndUnknowChrOVERLAP %>% head

## parents: 
uniteCovHALF_G1_woSexAndUnknowChrOVERLAP@treatment

## offspring:

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

