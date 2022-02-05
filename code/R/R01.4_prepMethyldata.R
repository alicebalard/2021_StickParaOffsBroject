## Alice Balard
## 4 Feb 2022
## NB: change the files to load each time R01.2 prep methylkit object is relaunched.
## Don't forget the timestamp to not mix up big files

## change path depending on the machine
if (machine=="apocrita"){
    mypath = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/"
} else if (machine=="mythinkpad"){
    mypath = "~/Documents/pro/Git/StickParaOffsBroject/gitignore/bigdata/"
}

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
base::load(paste0(mypath, "uniteCovALL_woSexAndUnknownChr4feb22.RData"))
## rename
fullMethylKitObj = uniteCovALL_woSexAndUnknowChr; rm(uniteCovALL_woSexAndUnknowChr)

## For further analyses: CpG covered in at least 2 individuals
## (Kostas took 2; Melanie is more stringent; let's see which makes sense)
base::load(paste0(mypath, "uniteCov2_woSexAndUnknownChr4feb22.RData"))
## rename
AL2MethylKitObj = uniteCov2_woSexAndUnknowChr ; rm(uniteCov2_woSexAndUnknowChr)

## If needed: CpG covered in at least 6 individuals
# base::load(paste0(mypath, "uniteCov6_woSexAndUnknownChr4feb22.RData"))

##### Create methylkit objects for different analyses:
# create a methylKit object with ONLY the parents - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_PAR=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_PAR$ID,
  treatment=fullMetadata_PAR$trtG1G2_NUM)
## methylBase object with 1772152 rows

## Remove positions with only NA
rm_NA_CpGs <- function(myunite){
    N=length(myunite@sample.ids)
    mydata=methylKit::getData(myunite)
    A=rowSums(is.na(mydata[grep("coverage", names(mydata))]))
    tab=table(A)
    tab[names(tab) %in% as.character(N)]#these CpGs have only NA
    ## CpGs to keep:
    mynonNAcols=which(rowSums(is.na(mydata[grep("coverage", names(mydata))]))!=N)
    ## SUBSET TO REMOVE FULL NA ROWS!!!!:
    myunite <- myunite[mynonNAcols]
    return(myunite)
}

uniteCov2_woSexAndUnknowChr_PAR <- rm_NA_CpGs(uniteCov2_woSexAndUnknowChr_PAR)

# create a methylKit object with ONLY the offspring - positions shared by ALL
uniteCovALL_woSexAndUnknowChr_OFF=reorganize(
  fullMethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)

## Remove positions with only NA
uniteCovALL_woSexAndUnknowChr_OFF <- rm_NA_CpGs(uniteCovALL_woSexAndUnknowChr_OFF)

# create a methylKit object with ONLY the offspring - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_OFF=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)
## methylBase object with 1772013 rows

## Remove positions with only NA
uniteCov2_woSexAndUnknowChr_OFF <- rm_NA_CpGs(uniteCov2_woSexAndUnknowChr_OFF)
