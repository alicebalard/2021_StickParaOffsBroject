#```{r focusFKBP3}
P = plotGeneTarget(myTargetGene = "FKBP3")
# Zoom in 1 
P1= P$plotGeneTarget + coord_cartesian(xlim=c(14053485,14053550)) + theme(legend.position = "none") +  ggtitle("") + labs(caption = "")
# Zoom in 2
P2=P$plotGeneTarget + coord_cartesian(xlim=c(14054790,14055010)) + theme(legend.position = "none") +  ggtitle("") + labs(caption = "")

# Plot gene and zooms
pdf("../../dataOut/FKBP3_DMS.pdf", width = 15, height = 5)
cowplot::plot_grid(cowplot::plot_grid(P1, P2, nrow = 1, ncol =2), P$plotGeneTarget, nrow = 2, ncol =1, rel_heights = c(1,2))
dev.off()

# Plot the different effects
DMSvec = P$myTargetGene_DMSdf$DMS[P$myTargetGene_DMSdf$effect %in% "G1"]

test = P$myTargetGene_DMSdf#[P$myTargetGene_DMSdf$effect %in% "G1",]


## Plot the most important CpG associated with significant axis 
# plotCpGexample_PCAsignif = function(DMSvec, mygene){
rm(test)

# Extract raw methylation values
raw=methylKit::select(uniteCov14_G2_woSexAndUnknowChrOVERLAP,
                      which(uniteCov14_G2_woSexAndUnknowChrOVERLAP$chr %in% test$chrom & 
                              uniteCov14_G2_woSexAndUnknowChrOVERLAP$start %in% sapply(strsplit(test$DMS," "), `[`, 2)))
dfPlot = data.frame(chr=raw$chr, pos=raw$end)
dfPlot = cbind(dfPlot, data.frame(percMethylation(raw)))
dfPlot = melt(dfPlot, id.vars = c("chr", "pos")) %>% dplyr::rename("SampleID" = "variable")
# Add sample group
dfPlot = merge(dfPlot, fullMetadata_OFFS[c("SampleID", "trtG1G2", "outcome", "patTrt", "brotherPairID")])
# Add effect
dfPlot$DMS = paste(dfPlot$chr, dfPlot$pos)
dfPlot = merge(dfPlot, test[c("effect", "DMS")])
# Plot
dfPlotSum = dfPlot %>% group_by(pos, outcome, patTrt, brotherPairID, effect) %>%
  dplyr::summarise(meanMeth=mean(value, na.rm=T)) %>% data.frame()

geom_rect(data = dfPlotSum, aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax =Inf, fill = effect), alpha = .05)
# make a df only for one value per BP to have a cool alpha
  

## All BP
ggplot()+
  # add effect:
  geom_rect(data = dfPlotSum, aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax =Inf, fill = effect), alpha = .05)+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))+#cb friendly palette
  # add points:
  geom_point(data = dfPlotSum, aes(x=outcome, y=meanMeth, col=patTrt, group=patTrt), size = 3) +
  geom_line(data = dfPlotSum, aes(x=outcome, y=meanMeth, col=patTrt, group=patTrt))+
  scale_color_manual("Paternal treatment", values = c("black", "red"))+
  # all data:
  geom_point(data = dfPlot, aes(x=outcome, y=value, col=patTrt, group=patTrt), alpha=.5) +
  facet_grid(brotherPairID~pos)+
  xlab("Offspring treatment")+
  ylab("Mean methylation")
