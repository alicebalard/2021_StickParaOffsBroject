# fication of
# ne by targeting methylated sites compressed between 1500bp
# stream from the transcription starting site (TSS) (Segonas et
#                                                    al., 2020). To identify nearest transcription starting site, GenomicRanges R package
# l., 2013) was used.
# 
#  enrichment and pathway analyses using the ENSEMBL
# ntify functional associations for both the top 2.5% highest
# he different parental treatments; and differentially hypermethylated
#  of over- representation of GO terms was evaluated using the
#  (Falcon & Gentleman, 2007). At a false discovery rate (FDR)
# 
#  functions were categorized as biological processes, molecular
# components. Visualization of gene functions was done with ggplot2
# 


############################## Identify Genes associated with positions 

############################## + GO terms --> TBC

# BiocManager::install("GOstats")
library("GOstats")
