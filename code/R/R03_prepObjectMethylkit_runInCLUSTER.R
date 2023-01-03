## MethylKit object preparation
## A. Balard
## 25th of August 2021 (updated Dec22)

# Each script sources the previous script of the pipeline if needed
# NB: comment when on Apocrita
source("R02_prepBSBOLTForMethylkit_runInCLUSTER.R")

## Load unitecov objects
base::load("../../gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovALL_woSexAndUnknowChr_20dec2022.RData") 
base::load("../../gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovALL_G1_woSexAndUnknowChr_20dec2022.RData") 
base::load("../../gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovALL_G2_woSexAndUnknowChr_20dec2022.RData") 
base::load("../../gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovHALF_G1_woSexAndUnknowChr_OVERLAPwG2_20dec2022.RData") 
base::load("../../gitignore/bigdata/05MethylKit/uniteCovObjects/uniteCovHALF_G2_woSexAndUnknowChr_OVERLAPwG1_20dec2022.RData")

rerun = FALSE

if (rerun == TRUE){

    ## Sources:
    ## https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html
    ## https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets
    ## /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

##### Load prepared dataset (in APOCRITA) #####
    dataPath="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/formatCG4methylKit"

    temp = list.files(path=dataPath,
                      pattern = ".CG4methylkit.txt",
                      full.names = T)

    ## Add metadata on treatments
    metadata <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/StickParaOffsBroject/data/raw\ data\ Joshka\ Kostas/Kostas_G2_info.xlsx")

    metadata$trtG1G2_NUM <- as.numeric(as.factor(metadata$trtG1G2))

### Make methylkit object
    myobj=methylKit::methRead(as.list(temp),
                              mincov=10,
                              sample.id=as.list(metadata$ID),
                              assembly="Gynogen_pchrom_assembly_all",
                              treatment=metadata$trtG1G2_NUM,
                              context="CpG")

#################################################
    ## We remove several samples from the raw dataset
    fullMetadata <- read.csv("../../data/fullMetadata127_Alice.csv") 

    ## create a new methylRawList object
    print("Remove unwanted samples")

    myobj=reorganize(
        myobj,
        sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID],
        treatment=metadata$trtG1G2_NUM[metadata$ID %in% fullMetadata$ID])

###############################
    ## Filtering and normalising ##
###############################

    ## Filtering based on coverage:
                                        # It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
    print("Filter")
    filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)

    ## normalise the coverage
    print("Normalise")
    normFil.myobj=normalizeCoverage(filtered.myobj)

#####################
    ## MERGING SAMPLES ##
#####################
    ##In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples.

    table(metadata$trtG1G2[metadata$ID %in% fullMetadata$ID])
                                        #   Control  E_control  E_exposed    Exposed NE_control NE_exposed 
                                        #        8         28         28         8         28         27 

    print("Add CpG present in ALL individuals")
    uniteCovALL= methylKit::unite(normFil.myobj, mc.cores=8)
    uniteCovALL=as(uniteCovALL,"methylBase")

    ## In all PARENTS:
    uniteCovALL_G1 = reorganize(
        normFil.myobj,
        sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "P"],
        treatment=metadata$trtG1G2_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "P"])

    uniteCovALL_G1 = methylKit::unite(uniteCovALL_G1, mc.cores=8) # try with 8 cores
    uniteCovALL_G1 = as(uniteCovALL_G1,"methylBase")

    ## In all OFFSPRING:
    uniteCovALL_G2 = reorganize(
        normFil.myobj,
        sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"],
        treatment=metadata$trtG1G2_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"])

    uniteCovALL_G2 = methylKit::unite(uniteCovALL_G2, mc.cores=8)# try with 8 cores
    uniteCovALL_G2 = as(uniteCovALL_G2,"methylBase")

    ## Keep methylated CpG sites observed in at least 50% individual fish after filtering and normalising: 4 for parents, 14 for offsprings
    ## PARENTS
    uniteCovHALF_G1 = reorganize(
        normFil.myobj,
        sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "P"],
        treatment=metadata$trtG1G2_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "P"])

    uniteCovHALF_G1 = methylKit::unite(uniteCovHALF_G1, min.per.group=4L, mc.cores=8)# try with 8 cores
    uniteCovHALF_G1 = as(uniteCovHALF_G1,"methylBase")

    ## OFFSPRING
    uniteCovHALF_G2 = reorganize(
        normFil.myobj,
        sample.ids=metadata$ID[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"],
        treatment=metadata$trtG1G2_NUM[metadata$ID %in% fullMetadata$ID & metadata$Generat %in% "O"])

    uniteCovHALF_G2 = methylKit::unite(uniteCovHALF_G2, min.per.group=14L, mc.cores=8)# try with 8 cores
    uniteCovHALF_G2 = as(uniteCovHALF_G2,"methylBase")

########################################################################################
    ## Remove reads from sex chromosome X ("Gy_chrXIX") and unmapped contigs ("Gy_chrUn") ##
########################################################################################
    print("nbr CpG on sex chromosome of unmapped:")
    nrow(uniteCovALL[uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),])

    print("Keep CpG apart from sex chromosome XIX and unmapped (comprise Y chr)")
    uniteCovALL_woSexAndUnknowChr=uniteCovALL[!uniteCovALL$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
    uniteCovALL_G1_woSexAndUnknowChr=uniteCovALL_G1[!uniteCovALL_G1$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
    uniteCovALL_G2_woSexAndUnknowChr=uniteCovALL_G2[!uniteCovALL_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

    uniteCovHALF_G1_woSexAndUnknowChr=uniteCovHALF_G1[!uniteCovHALF_G1$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]
    uniteCovHALF_G2_woSexAndUnknowChr=uniteCovHALF_G2[!uniteCovHALF_G2$chr %in% c("Gy_chrXIX", "Gy_chrUn"),]

    ## for positions covered in half fish groups, consider only methylKit object (uniteCov) with OVERLAPPING positions in Parents G1 and Offspring G2
    overlappingCpG_G1G2df <- findOverlaps(as(uniteCovHALF_G1_woSexAndUnknowChr,"GRanges"), 
                                          as(uniteCovHALF_G2_woSexAndUnknowChr,"GRanges"))
    overlappingCpG_G1G2df <- data.frame(overlappingCpG_G1G2df)

    uniteCovHALF_G1_woSexAndUnknowChrOVERLAP <- uniteCovHALF_G1_woSexAndUnknowChr[overlappingCpG_G1G2df$queryHits,]
    uniteCovHALF_G2_woSexAndUnknowChrOVERLAP <- uniteCovHALF_G2_woSexAndUnknowChr[overlappingCpG_G1G2df$subjectHits,]

    print("nbr CpG shared by all 127 samples:")
    length(uniteCovALL_woSexAndUnknowChr$chr)
    print("nbr CpG shared by all parents:")
    length(uniteCovALL_G1_woSexAndUnknowChr$chr)
    print("nbr CpG shared by all offsprings:")
    length(uniteCovALL_G2_woSexAndUnknowChr$chr)
    print("nbr CpG shared by 50% (4) of parents per trt group (overlapping with G2):")
    length(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP$chr)
    print("nbr CpG shared by 50% (14) of offsprings per trt group (overlapping with G1):")
    length(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP$chr)

######################################################################################################
    ## Save outcomes
    save(uniteCovALL_woSexAndUnknowChr,
         file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_woSexAndUnknowChr_20dec2022.RData")
    save(uniteCovALL_G1_woSexAndUnknowChr,
         file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G1_woSexAndUnknowChr_20dec2022.RData")
    save(uniteCovALL_G2_woSexAndUnknowChr,
         file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovALL_G2_woSexAndUnknowChr_20dec2022.RData")
    save(uniteCovHALF_G1_woSexAndUnknowChrOVERLAP,
         file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovHALF_G1_woSexAndUnknowChr_OVERLAPwG2_20dec2022.RData")
    save(uniteCovHALF_G2_woSexAndUnknowChrOVERLAP,
         file = "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/output/uniteCovHALF_G2_woSexAndUnknowChr_OVERLAPwG1_20dec2022.RData")

## NB: move outcome to gitignore to use on different machines
}

##################### Previous tests with ALL numbers of fish 1 to 12:
## we kept for downstream analyses all CpG sites present in at least 1 to 12 individuals per group, or in all individuals:
# print("Unite and store in a list")
# mylist_uniteCov=list()
# for (i in 6:12L){ # done for 1 to 5, then 6 to 12
#     uniteCov=unite(normFil.myobj, min.per.group=i, mc.cores=8)# try with 8 cores
#     uniteCov=as(uniteCov,"methylBase")
#     name=paste0("uniteCov_", as.character(i))
#     mylist_uniteCov[[name]]=uniteCov
# }

## Add CpG present in ALL individuals
# uniteCov=unite(normFil.myobj, mc.cores=8)
# uniteCov=as(uniteCov,"methylBase")
# mylist_uniteCov[["uniteCov_ALL"]]=uniteCov

# CpGALL=length(mylist_uniteCov$uniteCov_ALL$coverage1) # 47238

# Idea: plot number of retained CpG site by nbr of individuals sharing these CpG sites. Preparing file for that:
# print("Make DF")
# CpGDF1_5=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF1_5=t(CpGDF1_5)
# CpGDF1_5=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF1_5))),
#                     NbrCpG=CpGDF1_5[,1])
# 
# CpGDF6_12=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF6_12=t(CpGDF6_12)
# CpGDF6_12=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF6_12))),
#                      NbrCpG=CpGDF6_12[,1])
# 
# CpGDF=rbind(CpGDF1_5, CpGDF6_12)

# print("Save object for plotting")
# save(CpGDF, file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/CpGDF.RData")

# print("Save plot")
# plotCpGshared <- ggplot(CpGDF, aes(x=NbrIndMin, y=NbrCpG))+
#     geom_smooth(se = F, col = "red")+
#     geom_smooth(method = "lm", se = F, col = "black") +
#     geom_point() +
#     scale_x_continuous("Number of individual fish per treatment group sharing the same methylated CpG sites",
#                        labels = as.character(CpGDF$NbrIndMin), breaks = CpGDF$NbrIndMin)+
#     scale_y_continuous("Number of shared methylated CpG sites") +
#     theme_bw() +
#     geom_hline(yintercept=CpGALL)
# 
# plotCpGshared
# 
# pdf(file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/plotCpGshared.pdf")
# plotCpGshared
# dev.off()
