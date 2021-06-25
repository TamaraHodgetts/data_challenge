# Data challenge
# The directory is a folder called data
setwd("~/data_challenge/data (1)/data")

# import data
atacseq_design <- read.delim(file = "atacseq_design.txt", header = TRUE)
atacseq_peak_counts <- read.delim(file = "atacseq_peak_counts.txt", header = TRUE)
atacseq_peaks_bed <- read.delim(file = "atacseq_peaks.bed", header = FALSE)
rnaseq_annotation <- read.delim(file = "rnaseq_annotation.txt", header = TRUE)
rnaseq_design <- read.delim(file = "rnaseq_design.txt", header = TRUE)
rnaseq_gene_counts <- read.delim(file = "rnaseq_gene_counts.txt", header = TRUE)

# save the files in the proper format
save(atacseq_design,file="atacseq_design.Rda")
save(atacseq_peak_counts,file="atacseq_peak_counts.Rda")
save(atacseq_peaks_bed,file="atacseq_peaks_bed.Rda")
save(rnaseq_annotation,file="rnaseq_annotation.Rda")
save(rnaseq_design,file="rnaseq_design.Rda")
save(rnaseq_gene_counts,file="rnaseq_gene_counts.Rda")

# I'll probably have to normalise with respet to transcript length for RNA and ATAC-seq data 
summary(rnaseq_annotation$width)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8     527    2575   23267   17541 2304640 

# Reformat col 1 atacseq_peak_counts with start and end and chrom columns
install.packages("tidyverse")
library(dplyr)
library(tidyverse)
library(magrittr)
install.packages("fansi")
library(fansi)

library(tidyr)
install.packages("utf8")
library(utf8)
# ensure tidyr, utf8 and dplyr are installed and loaded

atacseq_peak_counts <- atacseq_peak_counts %>%
  separate(peakid, into = c("chr", "start-end"), sep = ":")
# separates the first column

atacseq_peak_counts <- atacseq_peak_counts %>%
  separate(`start-end`, into = c("start", "end"), sep = "-")
# separates the remaining column in 2

atacseq_peak_counts <- atacseq_peak_counts %>% 
  unite("peakid", start:end, sep = "-", remove = FALSE) %>%
  unite("peakid", chr:peakid, sep = ":", remove = FALSE)


save(atacseq_peak_counts,file="atacseq_peak_counts.Rda")

# Annotations
summary(rnaseq_annotation$gene_biotype)

levels(rnaseq_annotation$gene_biotype)
#[1] "3prime_overlapping_ncRNA"           "antisense"                          "bidirectional_promoter_lncRNA"      "IG_C_gene"                         
#[5] "IG_C_pseudogene"                    "IG_D_gene"                          "IG_J_gene"                          "IG_J_pseudogene"                   
#[9] "IG_pseudogene"                      "IG_V_gene"                          "IG_V_pseudogene"                    "lincRNA"                           
#[13] "macro_lncRNA"                       "miRNA"                              "misc_RNA"                           "Mt_rRNA"                           
#[17] "Mt_tRNA"                            "non_coding"                         "polymorphic_pseudogene"             "processed_pseudogene"              
#[21] "processed_transcript"               "protein_coding"                     "pseudogene"                         "ribozyme"                          
#[25] "rRNA"                               "scaRNA"                             "scRNA"                              "sense_intronic"                    
#[29] "sense_overlapping"                  "snoRNA"                             "snRNA"                              "sRNA"                              
#[33] "TEC"                                "TR_C_gene"                          "TR_D_gene"                          "TR_J_gene"                         
#[37] "TR_J_pseudogene"                    "TR_V_gene"                          "TR_V_pseudogene"                    "transcribed_processed_pseudogene"  
#[41] "transcribed_unitary_pseudogene"     "transcribed_unprocessed_pseudogene" "unitary_pseudogene"                 "unprocessed_pseudogene"            
#[45] "vaultRNA"
# Peaks have been annotated for names and function.

# Loading packages
install.packages('BiocManager')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")



BiocManager::install(c("ATACseqQC", "BiFET", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg19"))

# Analysing the data
# Peak annotation: 
source('https://bioconductor.org/biocLite.R')
biocLite(ChIPseeker)

install.packages(ChIPseeker)
library(ChIPpeakAnno)


version



