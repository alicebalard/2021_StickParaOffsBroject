file <- read.csv("../data/multiqc_report_rawReads_generalstats.csv") 
mean(file$M.Seqs)

# 95% confidence interval
qnorm(0.975)*sd(file$M.Seqs)/sqrt(nrow(file))
