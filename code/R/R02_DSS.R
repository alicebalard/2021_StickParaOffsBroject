## Alice BALARD
## July 2021
## QMUL
#### DSS help: https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html

## Install DSS (Dispersion shrinkage for sequencing data) package
# DSS is an R library performing differntial analysis for count-based sequencing data. It detectes differentially expressed genes (DEGs) from RNA-seq, and differentially methylated loci or regions (DML/DMRs) from bisulfite sequencing (BS-seq). The core of DSS is a new dispersion shrinkage method for estimating the dispersion parameter from Gamma-Poisson or Beta-Binomial distributions.

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DSS")

