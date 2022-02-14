## Alice Balard
## 4 Feb 2022
## NB: change the files to load each time R01.2 prep methylkit object is relaunched.
## Don't forget the timestamp to not mix up big files

## change path depending on the machine
if (machine=="apocrita"){
    mypath = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/"
} else if (machine=="mythinkpad"){
    mypath = "~/Documents/pro/Git/StickParaOffsBroject/gitignore/bigdata/05MethylKit/output/"
}

## Load previously united data (all 6 treatments)
## uniteCovALL: CpG covered in ALL individuals (has no NAs, useful for exploratory clustering analyses)
base::load(paste0(mypath, "uniteCovALL_woSexAndUnknownChr_10feb22.RData"))
## rename
#fullMethylKitObj = uniteCovALL_woSexAndUnknowChr; rm(uniteCovALL_woSexAndUnknowChr)

## Load previously united data for PARENTS
# CpG covered in ALL parents
base::load(paste0(mypath, "uniteCovALL_G1_woSexAndUnknownChr_10feb22.RData"))
## CpG covered in 50% of individuals (6)
base::load(paste0(mypath, "uniteCov6_G1_woSexAndUnknownChr_10feb22.RData"))

## Load previously united data for OFFSPRINGS
## CpG covered in 50% of individuals (14)
base::load(paste0(mypath, "uniteCov14_G2_woSexAndUnknownChr_10feb22.RData"))
# CpG covered in ALL offsprings
base::load(paste0(mypath, "uniteCovALL_G2_woSexAndUnknownChr_10feb22.RData"))

## How many CpGs are covered in each dataset?
nrow(uniteCovALL_woSexAndUnknowChr) # 55 530
nrow(uniteCovALL_G1_woSexAndUnknowChr) # 148 860
nrow(uniteCovALL_G2_woSexAndUnknowChr) # 78 384
nrow(uniteCov6_G1_woSexAndUnknowChr) # 1 188 179
nrow(uniteCov14_G2_woSexAndUnknowChr) # 1 050 222
