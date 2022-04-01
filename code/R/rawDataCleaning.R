## Cleaning raw data, associate a father to each offsrping clutch
## 1st of April 2022

## Load data
dataG1 <- read.csv("../../data/raw data Joshka Kostas/data-G1-Kostas_JK.csv")
dataAll <- read.csv("../../data/raw data Joshka Kostas/Kostas_G2_info.csv")
key <- read.csv("../../data/raw data Joshka Kostas/key_matchV48_clutchID_G2_Alice.csv")

## NOTES:
# Family in the G2 file is the family for Generat P (parent) and the father’s family for Generat O (offspring). 
# The family for Generat O is actually clutch_ID in data-G1-Kostas. The link between the two is "key"

as.character(dataG1$Idplate) %in% dataAll$matchV48
as.character(dataG1$Idplate) 
dataAll$matchV48

## NAs in dataG1 are in fact 0
dataG1$number.of.worms[is.na(dataG1$number.of.worms)] <- 0

## Merge parental data and check that there was no mistakes on crucial variables
cleanDF <- merge(data.frame(matchV48 = dataG1$Idplate,
                            FamilyOfFather = dataG1$FamilyOfFather,
                            clutch.ID = dataG1$clutch.ID,
                            HK.W = dataG1$HK.W/10,
                            No.Worms = dataG1$number.of.worms,
                            Wnettofin = dataG1$W.net,
                            Slfin = dataG1$SL),
                 dataAll, 
                 all = TRUE)

## Add offspring clutch ID matching with matchV48
cleanDF <- merge(data.frame(matchV48 = as.character(key$matchV48),
                            clutch.ID.G2 = key$clutch.ID),
                 cleanDF,
                 all=TRUE)

cleanDF$clutch.ID <- gsub("NA", "", paste0(cleanDF$clutch.ID, cleanDF$clutch.ID.G2))

## Rm clutchIDG2
cleanDF <- cleanDF[!names(cleanDF) %in% "clutch.ID.G2"]

## Add father ID
cleanDF <- merge(na.omit(data.frame(FamilyOfFather = cleanDF$FamilyOfFather,
                                    clutch.ID = cleanDF$clutch.ID,
                                    fatherID = cleanDF$ID)),
                 cleanDF,
                 all = TRUE)

## Remove fathers that were not sequenced
cleanDF <- cleanDF[!is.na(cleanDF$Generat),]

## Split offsping and fathers to add fatherID and FamilyOfFather to offspring
cleanDFG1 <- cleanDF[cleanDF$Generat %in% "P",]
cleanDFG2 <- cleanDF[cleanDF$Generat %in% "O",]

cleanDFG2 <- merge(data.frame(FamilyOfFather = cleanDFG1$FamilyOfFather,
                              clutch.ID = cleanDFG1$clutch.ID,
                              fatherID = cleanDFG1$fatherID),
                   cleanDFG2, by=c("clutch.ID"))

# rm redundant empty columns
cleanDFG2 <- cleanDFG2[!names(cleanDFG2) %in% c("FamilyOfFather.y", "fatherID.y")]
names(cleanDFG2)[names(cleanDFG2) %in% c("FamilyOfFather.x", "fatherID.x")] <- c("FamilyOfFather", "fatherID")

## Bind G1 and G2
cleanDF <- rbind(cleanDFG1, cleanDFG2)

## rm fatherID from fathers to avoid confusion
cleanDF$fatherID[cleanDF$Generat %in% "P"] <- NA

## Export cleaned dataset
write.csv(cleanDF, "../../data/cleanedRawData144fishG1G2.csv", quote = F, row.names = F)
