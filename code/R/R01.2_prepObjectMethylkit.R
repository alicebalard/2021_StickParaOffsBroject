## MethylKit analysis
## A. Balard
## 25th of August 2021

## R/4.0.2 to run
library(methylKit)
library(readxl)

## Sources:
## https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html

# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets

# /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

##### Load prepared dataset #####
dataPath="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/MethylationCalling/Methylation_calling_splitted/formatCG4methylKit"

temp = list.files(path=dataPath,
                   pattern = ".CG4methylkit.txt",
                   full.names = T)

## Add metadata on treatments
metadata <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/GIT_StickParaOffsBroject/data/Kostas_G2_info.xlsx")

metadata$trtG1G2_NUM <- as.numeric(as.factor(metadata$trtG1G2))

### Make methylkit object
myobj=methylKit::methRead(as.list(temp),
               mincov=10,
               sample.id=as.list(metadata$ID),
               assembly="Gynogen_pchrom_assembly_all",
               treatment=metadata$trtG1G2_NUM,
               context="CpG")

getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

## Filtering based on coverage:
# It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.

filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

## normalise the coverage
normFil.myobj=normalizeCoverage(filtered.myobj)

###### remove S12 that has a very low coverage and weird fastQC plots
## create a new methylRawList object
normFil.myobj143=reorganize(normFil.myobj,
                            sample.ids=metadata$ID[!metadata$ID %in% "S12"],
                            treatment=metadata$trtG1G2_NUM[!metadata$ID %in% "S12"])

## MERGING SAMPLES: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples.
table(metadata$trtG1G2)

##### Control  E_control  E_exposed    Exposed NE_control NE_exposed 
#        14         29         29         14         30         28 

## we kept for downstream analyses all CpG sites present in at least ten individuals per group:
uniteCov10=unite(normFil.myobj143, min.per.group=10L)

uniteCov10_mem <- as(uniteCov10,"methylBase")

save(uniteCov10_mem, file= "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCov10.RData")

## For correlation analyses (i.e. to test batch effects), we need no NA, so CpG present in ALL samples
uniteCovALL=unite(normFil.myobj143)

uniteCovALL_mem <- as(uniteCovALL,"methylBase")

save(uniteCovALL_mem, file= "/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/uniteCovALL.RData")
     
######### End of run
