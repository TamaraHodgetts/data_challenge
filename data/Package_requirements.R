# List of required packages

install.packages("knitr")

install.packages("rmdformats")

install.packages("dplyr")

install.packages("DT")

install.packages("tidyr")

install.packages("ggplot2")

install.packages("magrittr")

install.packages("devtools")

source('https://bioconductor.org/biocLite.R')

## Needed for mac and Linux only biocLite(Rsubread) ##

biocLite(Rsamtools)

biocLite(GenomicAlignments)

biocLite(TxDb.Hsapiens.UCSC.hg19.knownGene)

biocLite(soGGi)

biocLite(rtracklayer)

biocLite(ChIPQC)

biocLite(ChIPseeker)

biocLite(rGREAT)

biocLite(limma)

biocLite(DESeq2)

biocLite(tracktables)

biocLite(clusterProfiler)

biocLite(org.Mm.eg.db)

biocLite(MotifDb)

biocLite(Biostrings)

biocLite(BSgenome.Hsapiens.UCSC.hg19)

# # Finally we need development version of soGGi (named here 1.10.4) # not
# version on Bioconductor (1.10.0)
devtools::install_github('ThomasCarroll/soGGi')