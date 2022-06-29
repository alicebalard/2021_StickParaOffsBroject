## Load and prepare different types of metadata
## 2 feb 2022

#####################################################
## Load methylation metadata merged with Kostas files
fullMetadata <- read.csv("../../data/fullMetadata135_Alice.csv")

## Change one term: NbrMethylatedCpG is the one calculated by BSBolt IN TOTAL
names(fullMetadata)[names(fullMetadata) %in% "NbrMethylatedCpG"] <- "NbrMethylatedCpG_global_BSBolt"

# give a numerical value to treatment, for Methylkit
fullMetadata$trtG1G2_NUM <- as.numeric(as.factor(fullMetadata$trtG1G2))

## relevel treatments for graphs
fullMetadata$trtG1G2 <- factor(as.factor(fullMetadata$trtG1G2), levels = c("Control", "Exposed","NE_control", "NE_exposed", "E_control", "E_exposed"  ))

## family as factor for models
fullMetadata$Family <- as.factor(fullMetadata$Family)

# paternal exposure
fullMetadata$PAT="Exposed father group"
fullMetadata$PAT[fullMetadata$trtG1G2 %in% c("Control", "NE_control", "NE_exposed")]="Control father group"

## Add brother pairs
fullMetadata$brotherPairID <- sapply(str_split(fullMetadata$clutch.ID, "_"), "[", 2 )
## avoid confusion with numeric
fullMetadata$brotherPairID <- paste0("BP",fullMetadata$brotherPairID)

########################
## Parents only metadata
fullMetadata_PAR <- fullMetadata[fullMetadata$Generat %in% "P",]

##########################
## Offspring only metadata
fullMetadata_OFFS <- fullMetadata[fullMetadata$Generat %in% "O",]
fullMetadata_OFFS$trtG1G2 <- droplevels(fullMetadata_OFFS$trtG1G2)

## Create variable for offsping and parents separated
fullMetadata_OFFS$offsTrt <- "controlO"
fullMetadata_OFFS$offsTrt[fullMetadata_OFFS$Tr %in% c("TT", "CT")] <- "infectedO"
fullMetadata_OFFS$patTrt <- "controlP"
fullMetadata_OFFS$patTrt[fullMetadata_OFFS$Tr %in% c("TC", "TT")] <- "infectedP"

## Sanity check
table(fullMetadata_OFFS$offsTrt, fullMetadata_OFFS$trtG1G2)
table(fullMetadata_OFFS$patTrt, fullMetadata_OFFS$trtG1G2)

## REORDER metadata by sample ID
fullMetadata = fullMetadata[order(as.numeric(gsub("S", "", fullMetadata$SampleID))),]
fullMetadata_PAR = fullMetadata_PAR[order(as.numeric(gsub("S", "", fullMetadata_PAR$SampleID))),]
fullMetadata_OFFS = fullMetadata_OFFS[order(as.numeric(gsub("S", "", fullMetadata_OFFS$SampleID))),]

## Output cleaned data for Eri:
# write.csv(x = fullMetadata, file = "../../gitignore/fullMetadata4Eri.csv", quote = F, row.names = F)
