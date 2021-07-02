if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install(c("ATACseqQC", "BiFET", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg19"))
# load all packages to lib
library(ATACseqQC)
library(BiFET)
library(ChIPpeakAnno)
library(MotifDb)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)


#hg38 versions
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "phastCons100way.UCSC.hg38"))
# installed 1st, 2nd but not third.

# Install packages for visualisation
library("BiocManager")
library("Rsubread")
library("Rsamtools")
library("ggplot2")
install.packages("devtools")
library("usethis")
library("devtools")
library("magrittr")

