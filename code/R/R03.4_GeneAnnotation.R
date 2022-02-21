Furthermore, the identification of
methylated genes was done by targeting methylated sites compressed between 1500bp
downstream and 500bp upstream from the transcription starting site (TSS) (Segonas et
                                                                          al., 2020). To identify nearest transcription starting site, GenomicRanges R package
v.1.30.0 (Lawrence et al., 2013) was used.

We performed functional enrichment and pathway analyses using the ENSEMBL
reference genome to identify functional associations for both the top 2.5% highest
methylation ratios of the different parental treatments; and differentially hypermethylated
genes. The significance of over- representation of GO terms was evaluated using the
GOstats R package v2.52 (Falcon & Gentleman, 2007). At a false discovery rate (FDR)
6
threshold of 0.05, gene functions were categorized as biological processes, molecular
functions and cellular components. Visualization of gene functions was done with ggplot2
v3.3.2 (Wickham, 2016).