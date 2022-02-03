## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
#load("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknownChr.RData")
base::load("../../gitignore/output/uniteCovALL_woSexAndUnknownChr.RData")
## rename
fullMethylKitObj = uniteCovALL_woSexAndUnknowChr; rm(uniteCovALL_woSexAndUnknowChr)

## For further analyses: CpG covered in at least 2 individuals
## (Kostas took 2; Melanie is more stringent; let's see which makes sense)
load("../../gitignore/output/uniteCov2_woSexAndUnknownChr.RData")
## rename
AL2MethylKitObj = uniteCov2_woSexAndUnknowChr ; rm(uniteCov2_woSexAndUnknowChr)

## If needed: CpG covered in at least 6 individuals
#load("../../gitignore/output/uniteCov6_woSexAndUnknowChr.RData")

##### Create methylkit objects for different analyses:
# create a methylKit object with ONLY the parents - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_PAR=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_PAR$ID,
  treatment=fullMetadata_PAR$trtG1G2_NUM)

## Remove positions with only NA
# methylBase object with 1772152 rows
test <- methylKit::getData(uniteCov2_woSexAndUnknowChr_PAR)
A <- rowSums(is.na(test[grep("coverage", names(test))]))
table(A)#600226 columns have only NAs! Remove that
mynonNAcols <- which(rowSums(is.na(test[grep("coverage", names(test))]))!=24)
## SUBSET TO REMOVE FULL NA ROWS!!!!:
uniteCov2_woSexAndUnknowChr_PAR <- uniteCov2_woSexAndUnknowChr_PAR[mynonNAcols]
rm(test, A,mynonNAcols)

# create a methylKit object with ONLY the offspring - positions shared by ALL
uniteCovALL_woSexAndUnknowChr_OFF=reorganize(
  fullMethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)

# create a methylKit object with ONLY the offspring - positions shared by AT LEAST 2 INDIVIDUALS
uniteCov2_woSexAndUnknowChr_OFF=reorganize(
  AL2MethylKitObj,
  sample.ids=fullMetadata_OFFS$ID,
  treatment=fullMetadata_OFFS$trtG1G2_NUM)