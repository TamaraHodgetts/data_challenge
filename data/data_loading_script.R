# Data challenge
# The directory is a folder called data
setwd("~/data_challenge/data_(1)/data")

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

# separates the first column


atacseq_peak_counts <- atacseq_peak_counts %>%
  separate(`start-end`, into = c("start", "end"), sep = "-")

# separates the remaining column in 2

atacseq_peak_counts <- atacseq_peak_counts %>%
  unite("peakid", start:end, sep = "-", remove = FALSE) %>%
  unite("peakid", chr:peakid, sep = ":", remove = FALSE)


# save(atacseq_peak_counts,file="atacseq_peak_counts.Rda")

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
# Install all dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
# install Biocmanager
# update the old packages.

# load the package
library(BiocManager)
# and "ATACseqQC", "BiFET", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
#"BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg38.knownGene",
#"phastCons100way.UCSC.hg19"

# install needed packages
BiocManager::install(c("ATACseqQC", "BiFET", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg19"))


# Analysing the data
# Peak annotation:
#library(GenomeInfoDb)

#BiocManager::install("biocLite")

#biocLite(ChIPseeker)

#library(ChIPpeakAnno)
# BiocManager::install("GenomeInfoDbData")
# library(GenomeInfoDbData)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# BiocManager::install("GO.db")
# library(GO.db)
# BiocManager::install("DO.db")
# library(DO.db)
# BiocManager::install("ChIPseeker")
# library(ChIPseeker)
# 

# peakfile <- system.file("extdata", "atacseq_peak_counts.txt", package="ChIPseeker")
# ATAC_peaks_annotated <- annotatePeak(peakfile, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
# # peak should be GRanges object or a peak file...


# Differential expression analysis
BiocManager::install("Rsubread")
library(Rsubread)
install.packages("dplyr")
library(dplyr)

# All the peaks
# atacseq_peak_counts[ , 3:8]
# atacseq_peak_counts[ , 3:8] > 0
# summary(atacseq_peak_counts[ , 3:8] > 0)
# dim(atacseq_peak_counts[ , 3:8] > 0)
# # 173285      6
# discretised_atac_seq_peaks_df <- as.data.frame(atacseq_peak_counts[, 3:8] > 0)
# 
# atac_occurences_summary <- summary(discretised_atac_seq_peaks_df)
# atac_occurences_summary
# s84             s85             s86             s93             s94            s95         
# Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode:logical  
# FALSE:1710      FALSE:2807      FALSE:690       FALSE:41        FALSE:28        TRUE:173285   
# TRUE :171575    TRUE :170478    TRUE :172595    TRUE :173244    TRUE :173257 

# The maximal values and other summary statistics differ among the samples

# I want to keep only those peaks that are present in at least 3 replicates.
# What rows are inconsistent between replicates?
consistent_peak_atac_seq_subset <- subset(atacseq_peak_counts, (atacseq_peak_counts$s84 > 0 & atacseq_peak_counts$s85 > 0 &atacseq_peak_counts$s86 > 0) | (atacseq_peak_counts$s93 > 0 & atacseq_peak_counts$s94 > 0 &atacseq_peak_counts$s95 > 0))
dim(consistent_peak_atac_seq_subset)
# 173284      8

dim(atacseq_peak_counts)
# 173285      8
# There is one peak that is not consistently measured in either the control or treatment group. 
# this will be filtered out.

# Peak annotation
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# atacseq_peak_counts_anno <- annotatePeak(atacseq_peak_counts, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)


# Differential accessibility

atacseq_peak_counts <- atacseq_peak_counts %>%
  separate(peakid, into = c("chr", "start-end"), sep = ":")


gr2 <- GRanges(
  seqnames = atacseq_peak_counts$chr,
  ranges = atacseq_peak_counts$`start-end`)

gr2
save(gr2,file="Gene_ranges_atac_seq.Rda")
# Now I can use that for downstream analysis

# myPeaks <- GRanges(
#   seqnames = atacseq_peak_counts$chr,
#   ranges = atacseq_peak_counts$`start-end`,
#   control_1 = atacseq_peak_counts$s84,
#   control_2 = atacseq_peak_counts$s85,
#   control_3 = atacseq_peak_counts$s86,
#   treated_1 = atacseq_peak_counts$s93,
#   treated_2 = atacseq_peak_counts$s94,
#   treated_3 = atacseq_peak_counts$s95,)

dim(consistent_peak_counts_atac_seq_subset)
# 173284
chr_vec <- rep("chr", 173284)
length(chr_vec)

consistent_peak_atac_seq_subset <- consistent_peak_atac_seq_subset %>%
  mutate(chr_vec = chr_vec) %>%
  unite("seqnames", chr_vec, chr, sep = "")

consistent_peak_atac_seq_subset

library(GenomicRanges)

myPeaks <- GRanges(
  seqnames = consistent_peak_atac_seq_subset$seqnames,
  ranges = consistent_peak_atac_seq_subset$`start-end`,
  control_1 = consistent_peak_atac_seq_subset$s84,
  control_2 = consistent_peak_atac_seq_subset$s85,
  control_3 = consistent_peak_atac_seq_subset$s86,
  treated_1 = consistent_peak_atac_seq_subset$s93,
  treated_2 = consistent_peak_atac_seq_subset$s94,
  treated_3 = consistent_peak_atac_seq_subset$s95,)

myPeaks

Group <- factor(c("control", "control", "control", "treated", "treated", "treated"))

# venn diagram overlap graphs
library(limma)

# for atac seq controls
venn_atac_seq_controls <- as.data.frame(elementMetadata(myPeaks)) %>% dplyr::select(starts_with("control")) %>% 
  vennDiagram(main = "Overlap for control group open regions")

pdf("Venn_diagram_atac_seq_control.pdf")
as.data.frame(elementMetadata(myPeaks)) %>% dplyr::select(starts_with("control")) %>% 
  vennDiagram(main = "Overlap for control group open regions")
dev.off()

# for atac seq treated
venn_atac_seq_treated <- as.data.frame(elementMetadata(myPeaks)) %>% dplyr::select(starts_with("treated")) %>% 
  vennDiagram(main = "Overlap for treated group open regions")

pdf("Venn_diagram_atac_seq_treated.pdf")
as.data.frame(elementMetadata(myPeaks)) %>% dplyr::select(starts_with("treated")) %>% 
  vennDiagram(main = "Overlap for treated group open regions")
dev.off()

# There is higher variability in the chromatin accessibility profiles in the control groups as compared to the treated groups

# Subset those peaks that are present in at least 3 replicates belonging to the same class (i.e. control or treated)
# Filter away from the myPeaks object those peaks that are not consistently present in either the
# 3 control or the 3 treatment samples



# library(Rsubread)
# 
# occurrences <- elementMetadata(myPeaks) %>% as.data.frame %>%  rowSums
# occurrences
# dim(occurences)
# # 173285      6
# 
# ATAC_seq_occurences_df <- as.data.frame(occurences)
# # This is the same as the peak counts dataframe
# table(occurrences) %>% rev %>% cumsum
# 
# atac_seq_peaks_df
# rowSums(atac_seq_peaks_df[, 6:11])
# 
# table(occurrences) %>% rev %>% cumsum
# 
# rm(ATAC_seq_occurences_df)
# rm(occurences)

# PCA for overlaps
library(tidyr)

myPeaks

myPlot <- as.data.frame(elementMetadata(myPeaks)) %>% as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_//d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5)

myPlot

pdf("PCA of atac-seq data (treatment and control).pdf")
as.data.frame(elementMetadata(myPeaks)) %>% as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_//d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5)
dev.off()


# Differential accessibility analysis

atac_seq_peaks_df <- as.data.frame(myPeaks)
view(atac_seq_peaks_df)

dim(atac_seq_peaks_df)
# 173284     11

dim(consistent_peak_counts_atac_seq_subset)
#  173284     7

# Number of rows is the same as the atacseq_peak_counts object


# Deseq2 analysis
library(DESeq2)
metaData <- data.frame(Group, row.names = colnames(myPeaks))
metaData


dim(metaData)
# 6  1
dim(atac_seq_peaks_df)
# 173284     11
dim(consistent_peak_counts_atac_seq_subset)
# 173284      8

subsetted_atac_seq_peak_counts <- consistent_peak_atac_seq_subset[, 3:8]
subsetted_atac_seq_peak_counts
dim(subsetted_atac_seq_peak_counts)


atacDDS <- DESeqDataSetFromMatrix(subsetted_atac_seq_peak_counts, metaData, ~Group, rowRanges = myPeaks)
atacDDS
# To use the function DESeqDataSetFromMatrix, ncol(countData) must equal nrow(colData) 

# DESeqDataSetFromMatrix(
#   countData,
#   colData,
#   design,
#   tidy = FALSE,
#   ignoreRank = FALSE,
#   ...
# )

atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)
plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))

plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog)) + 
  labs(
    title = paste(
      "PCA plot of ATAC-seq dataset"))


# save the file 
pdf("ATAC-seq_PCA.pdf")
plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog)) + 
  labs(
    title = paste(
      "PCA plot of ATAC-seq dataset"))
dev.off()


# Differential open analysis
library(DESeq2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tracktables)


atac_seq_control_treatment <- results(atacDDS, c("Group", "control", "treated"), format = "GRanges")
atac_seq_control_treatment <- atac_seq_control_treatment[order(atac_seq_control_treatment$pvalue)]


atac_seq_control_treatment

save(atac_seq_control_treatment, file="atac_seq_deseq2.Rda")


# Subsetting those regions of open chromatin within promoters
# and creating a table to be able to view the results in IGV
# using the makebedtable function from the tracktables package
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# toOverLap <- promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, 500, 500)
# toOverLap

toOverLap_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "gene_id", filter=NULL, single.strand.genes.only=TRUE)
toOverLap_genes


# subsetting regions of significantly enriched open chromatin
 atac_seq_control_treatment_significant <- atac_seq_control_treatment[(!is.na(atac_seq_control_treatment$padj) & 
                                                                         atac_seq_control_treatment$padj < 0.01) & atac_seq_control_treatment %over% toOverLap_genes, ]

 atac_seq_control_treatment_significant
 
# 11015
# The filtering is the same, regartless of whether is uses the genes or promoters for the gene annotation.  
 
# Warning message:
 # In .Seqinfo.mergexy(x, y) :
 #   Each of the 2 combined objects has sequence levels not in the other:
 #   - in 'x': chrGL000008.2, chrGL000009.2, chrGL000194.1, chrGL000195.1, chrGL000205.2, chrGL000208.1, chrGL000213.1, chrGL000214.1, chrGL000216.2, chrGL000218.1, chrGL000219.1, chrGL000220.1, chrGL000221.1, chrGL000224.1, chrGL000225.1, chrGL000226.1, chrKI270303.1, chrKI270304.1, chrKI270310.1, chrKI270311.1, chrKI270312.1, chrKI270317.1, chrKI270320.1, chrKI270330.1, chrKI270333.1, chrKI270334.1, chrKI270336.1, chrKI270337.1, chrKI270340.1, chrKI270384.1, chrKI270391.1, chrKI270395.1, chrKI270396.1, chrKI270411.1, chrKI270412.1, chrKI270418.1, chrKI270422.1, chrKI270435.1, chrKI270438.1, chrKI270442.1, chrKI270465.1, chrKI270466.1, chrKI270467.1, chrKI270512.1, chrKI270517.1, chrKI270519.1, chrKI270538.1, chrKI270544.1, chrKI270580.1, chrKI270581.1, chrKI270584.1, chrKI270587.1, chrKI270589.1, chrKI270590.1, chrKI270591.1, chrKI270706.1, chrKI270709.1, chrKI270711.1, chrKI270712.1, chrKI270713.1, chrKI270714.1, chrK [... truncated]                                                              suppressWarnings() to suppress this warning.)

makebedtable(atac_seq_control_treatment_significant, "atac_seq_control_treatment_significant.html", basedirectory = getwd())


# I can subset the statistically significant ATAC-seq peaks

# atac_seq_control_treatment_significant <- atac_seq_control_treatment[(!is.na(atac_seq_control_treatment$padj) & 
#                                                             atac_seq_control_treatment$padj < 0.01), ]
# 
# atac_seq_control_treatment_significant

# and save the deseq results.
save(atac_seq_control_treatment_significant,file="significant_deseq2_atac-seq_peaks.Rda")


# The same can be done with the RNA-seq data


# myPeaks_rna_seq_all <- GRanges(
#   seqnames = rnaseq_gene_counts$featureid,
#   ranges = rnaseq_annotation_united$`start-end`,
#   control_1 = rnaseq_gene_counts$s69,
#   control_2 = rnaseq_gene_counts$s70,
#   control_3 = rnaseq_gene_counts$s71,
#   treated_1 = rnaseq_gene_counts$s75,
#   treated_2 = rnaseq_gene_counts$s76,
#   treated_3 = rnaseq_gene_counts$s77,)
# 
# myPeaks_rna_seq_all



# venn diagram overlap graphs
library(limma)

# for rna seq controls
venn_rna_seq_controls <- as.data.frame(elementMetadata(myPeaks_rna_seq_all)) %>% dplyr::select(starts_with("control")) %>% 
  vennDiagram(main = "Overlap for control group transcribed regions")

pdf("Venn_diagram_rna_seq_control.pdf")
as.data.frame(elementMetadata(myPeaks_rna_seq_all)) %>% dplyr::select(starts_with("control")) %>% 
  vennDiagram(main = "Overlap for control group transcribed regions")
dev.off()

# for rna seq treated
venn_rna_seq_treated <- as.data.frame(elementMetadata(myPeaks_rna_seq_all)) %>% dplyr::select(starts_with("treated")) %>% 
  vennDiagram(main = "Overlap for treated group transcribed regions")

pdf("Venn_diagram_rna_seq_treated.pdf")
as.data.frame(elementMetadata(myPeaks_rna_seq_all)) %>% dplyr::select(starts_with("treated")) %>% 
  vennDiagram(main = "Overlap for treated group transcribed regions")
dev.off()

# You can use these kinds of graphs for IDR analysis

# dim(consistent_peak_counts_atac_seq_subset)
# # 173284
# chr_vec <- rep("chr", 173284)
# length(chr_vec)
# 
# consistent_peak_counts_atac_seq_subset <- consistent_peak_counts_atac_seq_subset %>%
#   mutate(chr_vec = chr_vec) %>%
#   unite("seqnames", chr_vec, chr, sep = "")
# 
# consistent_peak_counts_atac_seq_subset
# 
# myPeaks <- GRanges(
#   seqnames = consistent_peak_counts_atac_seq_subset$seqnames,
#   ranges = consistent_peak_counts_atac_seq_subset$`start-end`,
#   control_1 = consistent_peak_counts_atac_seq_subset$s84,
#   control_2 = consistent_peak_counts_atac_seq_subset$s85,
#   control_3 = consistent_peak_counts_atac_seq_subset$s86,
#   treated_1 = consistent_peak_counts_atac_seq_subset$s93,
#   treated_2 = consistent_peak_counts_atac_seq_subset$s94,
#   treated_3 = consistent_peak_counts_atac_seq_subset$s95,)
# 
# myPeaks
# 
# Group <- factor(c("control", "control", "control", "treated", "treated", "treated"))
# 
# rnaseq_gene_counts
# 
# rnaseq_annotation_united$chr
install.packages("dplyr")
library(dplyr)

rnaseq_annotation %>%
  unite("start-end", start:end, sep = "-", remove = FALSE)

rnaseq_annotation_united <- rnaseq_annotation %>%
  unite("start-end", start:end, sep = "-", remove = FALSE)

dim(rnaseq_annotation_united)

chr_vec_rna <- rep("chr", 58051)

rnaseq_annotation_united <- rnaseq_annotation_united %>%
  mutate(chr_vec_rna = chr_vec_rna) %>%
  unite("seqnames", chr_vec_rna, chr, sep = "")
rnaseq_annotation_united

# Checking that the nr of rows are the same
dim(rnaseq_annotation_united)
# 58051    33

dim(rnaseq_gene_counts)
# 58051     7

summary(rnaseq_annotation_united$featureid == rnaseq_gene_counts$featureid)
# Mode    TRUE 
# logical   58051
# The feature ids are identical

# Subsetting consistent RNA-seq peaks
consistent_peak_annotations_rna_seq_subset <- subset(rnaseq_annotation_united, (rnaseq_gene_counts$s69 > 0 & rnaseq_gene_counts$s70 > 0 &rnaseq_gene_counts$s71 > 0) | (rnaseq_gene_counts$s75 > 0 & rnaseq_gene_counts$s76 > 0 &rnaseq_gene_counts$s77 > 0))
dim(consistent_peak_annotations_rna_seq_subset)
#  20289    33
consistent_peak_annotations_rna_seq_subset$seqnames

consistent_peak_counts_rna_seq_subset <- subset(rnaseq_gene_counts, (rnaseq_gene_counts$s69 > 0 & rnaseq_gene_counts$s70 > 0 &rnaseq_gene_counts$s71 > 0) | (rnaseq_gene_counts$s75 > 0 & rnaseq_gene_counts$s76 > 0 &rnaseq_gene_counts$s77 > 0))
dim(consistent_peak_counts_rna_seq_subset)
#  20289    7

consistent_peak_annotations_rna_seq_subset
consistent_peak_counts_rna_seq_subset

summary(consistent_peak_annotations_rna_seq_subset$featureid == consistent_peak_counts_rna_seq_subset$featureid)
# Mode    TRUE 
# logical   20289

# There is one peak that is not consistently measured in either the control or treatment group. 
# this will be filtered out.


myPeaks_rna_seq <- GRanges(
  seqnames = consistent_peak_annotations_rna_seq_subset$seqnames,
  ranges = consistent_peak_annotations_rna_seq_subset$`start-end`,
  control_1 = consistent_peak_counts_rna_seq_subset$s69,
  control_2 = consistent_peak_counts_rna_seq_subset$s70,
  control_3 = consistent_peak_counts_rna_seq_subset$s71,
  treated_1 = consistent_peak_counts_rna_seq_subset$s75,
  treated_2 = consistent_peak_counts_rna_seq_subset$s76,
  treated_3 = consistent_peak_counts_rna_seq_subset$s77,)

myPeaks_rna_seq

Group <- factor(c("control", "control", "control", "treated", "treated", "treated"))


# Subsetting consistent peaks
# In the RNA-seq datasets, there is much more discordance between the sample replicates, as compared to the ATAC-seq data.
# I need to subset those peaks that are consistent among the control samples, and those that are consistent among the treated samples.
# Then, I can compare these consistent peaks between the control and treated samples

# PCA for overlaps
library(tidyr)

myPeaks_rna_seq

myPlot_rna_seq <- as.data.frame(elementMetadata(myPeaks_rna_seq)) %>% as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_//d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5)

myPlot_rna_seq

pdf("PCA of rna-seq data (treatment and control).pdf")
as.data.frame(elementMetadata(myPeaks_rna_seq)) %>% as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_//d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5)
dev.off()

# Differential accessibility analysis for RNA seq

# Deseq2 analysis
library(DESeq2)
metaData_rna_seq <- data.frame(Group, row.names = colnames(myPeaks_rna_seq))
metaData_rna_seq


dim(metaData_rna_seq)
# 6  1

dim(consistent_peak_counts_rna_seq_subset)
# 20289    7
dim(consistent_peak_annotations_rna_seq_subset)
# 20289     33

subsetted_rnaseq_gene_counts <- consistent_peak_counts_rna_seq_subset[, 2:7]
subsetted_rnaseq_gene_counts
dim(subsetted_rnaseq_gene_counts)
# 20289     6

rnaDDS <- DESeqDataSetFromMatrix(subsetted_rnaseq_gene_counts, metaData_rna_seq, ~Group, rowRanges = myPeaks_rna_seq)
rnaDDS


rnaDDS <- DESeq(rnaDDS)
rnaDDS

# class: DESeqDataSet 
# dim: 20289 6 
# metadata(1): version
# assays(4): counts mu H cooks
# rownames(20289): 1 3 ... 58041 58048
# rowData names(28): control_1 control_2 ... deviance maxCooks
# colnames(6): s69 s70 ... s76 s77
# colData names(2): Group sizeFactor

rna_Rlog <- rlog(rnaDDS)

plotPCA(rna_Rlog, intgroup = "Group", ntop = nrow(rna_Rlog)) + 
  labs(
    title = paste(
      "PCA plot of RNA-seq dataset"))
    

# save that file and the ATAC-seq equivalent from yesterday
pdf("RNA-seq_PCA.pdf")
plotPCA(rna_Rlog, intgroup = "Group", ntop = nrow(rna_Rlog)) + 
  labs(
    title = paste(
      "PCA plot of RNA-seq dataset"))
dev.off()

# Differential expression analysis
library(DESeq2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tracktables)

# rna_seq_control_treatment_1 <- results(rnaDDS, c("Group", "control", "treated"), format = "GRanges")
# rna_seq_control_treatment_2 <- results(rnaDDS, c("Group", "treated", "control"), format = "GRanges")
# rm(rna_seq_control_treatment_1, rna_seq_control_treatment_2)

rna_seq_control_treatment <- results(rnaDDS, c("Group", "control", "treated"), format = "GRanges")

rna_seq_control_treatment <- rna_seq_control_treatment[order(rna_seq_control_treatment$pvalue)]
rna_seq_control_treatment

summary(rna_seq_control_treatment$padj)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.0000  0.0030  0.2454  0.3532  0.6854  0.9999    2362  

# What causes these NA values?

summary(atac_seq_control_treatment$padj)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1327  0.4267  0.4399  0.7257  1.0000 

save(rna_seq_control_treatment, file="rna_seq_deseq2.Rda")


# I can subset the statistically significant RNA-seq peaks

# Subsetting those regions of expressed chromatin within promoters
# and creating a table to be able to view the results in IGV
# using the makebedtable function from the tracktables package
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#toOverLap_promoters <- promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, 500, 500)
#toOverLap_promoters

toOverLap_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "gene_id", filter=NULL, single.strand.genes.only=TRUE)
toOverLap_genes


# subsetting regions of significantly enriched open chromatin (in a 1000 bp distance from a promotor)
rna_seq_control_treatment_significant <- rna_seq_control_treatment[(!is.na(rna_seq_control_treatment$padj) & 
                                                                      rna_seq_control_treatment$padj < 0.01) & rna_seq_control_treatment %over% toOverLap_genes, ]

rna_seq_control_treatment_significant
# 5077
makebedtable(rna_seq_control_treatment_significant, "rna_seq_control_treatment_significant.html", basedirectory = getwd())

save(rna_seq_control_treatment_significant,file="significant_deseq2_rna-seq_peaks.Rda")


# Annotation for differential atac-seq peaks
# Gene annotation
install.packages(cluserProfiler)
BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(ChIPseeker)

rna_seq_control_treatment_significant
atac_seq_control_treatment_significant

?annotatePeak

anno_atac_peaks <- annotatePeak(atac_seq_control_treatment_significant, tssRegion=c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_atac_peaks
# Annotated peaks generated by ChIPseeker
#11015/11015  peaks were annotated
# Genomic Annotation Summary:
#   Feature  Frequency
# 8   Promoter (<=1kb) 12.5646845
# 9   Promoter (1-2kb)  5.4380390
# 10  Promoter (2-3kb)  5.2473899
# 4             5' UTR  0.4993191
# 3             3' UTR  3.7948252
# 1           1st Exon  1.2074444
# 6         Other Exon  5.0385837
# 2         1st Intron 19.5097594
# 7       Other Intron 45.3835679
# 5  Distal Intergenic  1.3163867

rna_seq_control_treatment_significant_df

anno_rna_peaks <- annotatePeak(rna_seq_control_treatment_significant, tssRegion=c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_rna_peaks
# Annotated peaks generated by ChIPseeker
# 5077/5077  peaks were annotated
# Genomic Annotation Summary:
#   Feature   Frequency
# 5 Promoter (<=1kb) 98.30608627
# 6 Promoter (1-2kb)  0.25605673
# 7 Promoter (2-3kb)  0.25605673
# 3           5' UTR  0.15757337
# 2           3' UTR  0.53181012
# 1         1st Exon  0.39393343
# 4       Other Exon  0.09848336

# If I reduce the range, nearly all is promoter associated. 
# It's and Chip-seq function, applied to ATAC- and RNA-seq data (with generally shorter reads)

rna_seq_control_treatment_significant_df <- as.data.frame(rna_seq_control_treatment_significant)
rna_seq_control_treatment_significant_df

atac_seq_control_treatment_significant_df <- as.data.frame(atac_seq_control_treatment_significant)
atac_seq_control_treatment_significant_df

rna_seq_control_treatment_significant_df$start[1]
# 118607643

?results
# results extracts a result table from a DESeq analysis giving base means across samples, log2 fold 
# changes, standard errors, test statistics, p-values and adjusted p-values; resultsNames returns the 
# names of the estimated effects (coefficents) of the model; removeResults returns a DESeqDataSet object 
# with results columns removed.

# Visualising the gene annotation

plotAnnoPie(anno_atac_peaks)

pdf("Pie_chart_atac_seq_peak_annotation.pdf")
plotAnnoPie(anno_atac_peaks)
dev.off()

plotAnnoPie(anno_rna_peaks)

pdf("Pie_chart_rna_seq_peak_annotation.pdf")
plotAnnoPie(anno_rna_peaks)
dev.off()


plotAnnoBar(anno_rna_peaks)
plotAnnoBar(anno_atac_peaks)

install.packages("ggupset")
library(ggupset)

upsetplot(anno_rna_peaks)
upsetplot(anno_atac_peaks)

?upsetplot
# There's another upsetplot in the enrichplot library

pdf("Upset_plot_atac_seq_peak_annotation.pdf")
upsetplot(anno_atac_peaks)
dev.off()

pdf("Upset_plot_rna_seq_peak_annotation.pdf")
upsetplot(anno_rna_peaks)
dev.off()

# Subsetting those RNA-seq peaks that are transcriptionally enriched
# 
# atac_seq_control_treatment_significant
# rna_seq_control_treatment_significant
# 
# # Use a filtering join
# rna_seq_control_treatment_significant
# 
# rna_seq_control_treatment_significant_df
# dim(rna_seq_control_treatment_significant_df)
# # Uniting the start and end
# rna_seq_control_treatment_significant_df <- rna_seq_control_treatment_significant_df %>%
#   unite("start-end", start:end, sep = "-", remove = FALSE) 
# 
# rna_seq_control_treatment_significant_df
# view(rna_seq_control_treatment_significant_df)
# 
# rnaseq_annotation_united
# view(rnaseq_annotation_united)
# 
# # Joining the dfs
# joined_df <- semi_join(rnaseq_annotation_united, rna_seq_control_treatment_significant_df, by = "seqnames", "start", "end")
# joined_df
# 
# dim(rnaseq_annotation_united)
# # 58051    33
# 
# dim(rna_seq_control_treatment_significant_df)
# # 6363   11
# 
# dim(joined_df)
# # 57955    33
# 
# # TSSs
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), fix = "start", 1)
# TSSs


# Visualising distance or reads to the TSS
pdf("TSS_distance_RNA_peaks.pdf")
plotDistToTSS(anno_rna_peaks,
              title="Distribution of RNA-peaks/nrelative to TSS")
dev.off()

pdf("TSS_distance_ATAC_peaks.pdf")
plotDistToTSS(anno_atac_peaks,
              title="Distribution of ATAC-peaks/nrelative to TSS")
dev.off()

summary(anno_rna_peaks_df$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 65    8946   23949   46699   57223 1193736 

summary(consistent_peak_annotations_rna_seq_subset$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 50    3706   14698   40766   43826 2304640

summary(rnaseq_annotation$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8     527    2575   23267   17541 2304640 

# ANalysing differential accessibility expression analysis results and gene annotation
anno_atac_peaks_df <- as.data.frame(anno_atac_peaks)
view(anno_atac_peaks_df)

anno_rna_peaks_df <- as.data.frame(anno_rna_peaks)
view(anno_rna_peaks_df)

# Joining the dfs
summary(colnames(anno_atac_peaks_df) == colnames(anno_rna_peaks_df))
# Mode    TRUE 
# logical      20 

colnames(anno_atac_peaks_df)
# [1] "seqnames"       "start"          "end"            "width"          "strand"         "baseMean"       "log2FoldChange" "lfcSE"         
# [9] "stat"           "pvalue"         "padj"           "annotation"     "geneChr"        "geneStart"      "geneEnd"        "geneLength"    
# [17] "geneStrand"     "geneId"         "transcriptId"   "distanceToTSS" 

joined_df_transcriptId <- semi_join(anno_rna_peaks_df, anno_atac_peaks_df, by = "transcriptId")
joined_df_transcriptId
dim(joined_df_transcriptId)
# 446  20

# joined_df_gene <- semi_join(anno_rna_peaks_df, anno_atac_peaks_df, by = c("geneStart", "geneEnd"))
# joined_df_gene
# dim(joined_df_gene)
# 768  20

joined_df_geneId <- semi_join(anno_rna_peaks_df, anno_atac_peaks_df, by = c("geneId", "transcriptId"))
joined_df_geneId
dim(joined_df_geneId)
# 446  20

save(joined_df_geneId,file="Joined_ATAC_and_RNA_seq_df.Rda")

# 

significant_transcript_Ids <- joined_df_geneId$transcriptId
significant_gene_Ids <- joined_df_geneId$geneId

joined_df_geneId

# Save this data to a file and search it in eDirect

# Q2: Identifying enriched transcriptional programs
BiocManager::install("ReactomePA")
library(ReactomePA)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

dim(anno_rna_peaks_df)
#5077   20
dim(joined_df_geneId)
# 446  20
summary(colnames(anno_rna_peaks_df) == colnames(joined_df_geneId))
# Mode    TRUE 
# logical   20

rm(pathway1)
pathways_atac <- enrichPathway(as.data.frame(anno_atac_peaks)$geneId)
pathways_atac
head(pathways_atac, 2)
#      ID       Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# R-HSA-112316   R-HSA-112316   Neuronal System  166/2768 410/10856 8.954224e-12 1.347611e-08 1.144256e-08
# R-HSA-9013149 R-HSA-9013149 RAC1 GTPase cycle   85/2768 185/10856 1.066587e-09 8.026070e-07 6.814932e-07

pathways_atac_df <- as.data.frame(pathways_atac)

pathways_rna <- enrichPathway(as.data.frame(anno_rna_peaks)$geneId)
pathways_rna
head(pathways_rna, 2)
# ID               Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# R-HSA-187037 R-HSA-187037 Signaling by NTRK1 (TRKA)   61/2489 115/10856 1.837212e-12 2.750306e-09 2.084752e-09
# R-HSA-166520 R-HSA-166520        Signaling by NTRKs   64/2489 134/10856 2.006324e-10 1.501734e-07 1.138325e-07

pathways_rna_df <- as.data.frame(pathways_rna)
pathways_atac_df$ID


pdf("RNA_seq_pathway_enrichment.pdf")
dotplot(pathways_rna) +
  labs(
    title = paste(
      "Pathway enrichment analysis of significant RNA-seq peaks"))
dev.off()

pdf("RNA_seq_pathway_enrichment.pdf")
dotplot(pathways_atac)
labs(
  title = paste(
    "Pathway enrichment analysis of significant ATAC-seq peaks"))
dev.off()

# Identifying pathways in both programs
joined_pathways <- semi_join(pathways_rna_df, pathways_atac_df, by = "ID")

dim(joined_pathways)
# 38  9

# Saving the objects
save(joined_pathways, file="enriched_transcriptional_programs.Rda")
save(pathways_rna_df, file="rna_seq_pathway_enrichment.Rda")
save(pathways_atac_df, file="atac_seq_pathway_enrichment.Rda")


# Attempting to analyse GO term enrichment on the joined DF


joined_df_geneId

joined_df_geneId_united <- joined_df_geneId %>%
  unite("peakid", start:end, sep = "-", remove = FALSE)

myRanges <- GRanges(
  seqnames = joined_df_geneId_united$seqnames,
  ranges = joined_df_geneId_united$peakid,)

myRanges

anno_ranges <- annotatePeak(myRanges, tssRegion=c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_ranges
# Annotated peaks generated by ChIPseeker
# 446/446  peaks were annotated
# Genomic Annotation Summary:
#   Feature  Frequency
# 5 Promoter (<=1kb) 94.1704036
# 6 Promoter (2-3kb)  0.4484305
# 3           5' UTR  0.2242152
# 2           3' UTR  2.9147982
# 1         1st Exon  1.7937220
# 4       Other Exon  0.4484305


pdf("Pie_chart_rna_seq_peak_annotation.pdf")
plotAnnoPie(anno_ranges)
dev.off()

upsetplot(anno_ranges)

# Attempting to annotate the object
# anno_ranges_df <- as.data.frame(anno_ranges)
# 
# pathways_join <- enrichPathway(as.data.frame(anno_ranges)$geneId)
# pathways_join
# head(pathways_join, 2)


# joined_pathways
# joined_df_geneId
# 
# enrichPathway(joined_df_geneId$geneId)
# 
# GR_joined_pathways <- joined_df_geneId
# GR_joined_pathways$geneId
# 
# joined_pathways_enrichment <- enrichPathway(GR_joined_pathways$geneId)
# joined_pathways_enrichment_df <- as.data.frame(joined_pathways_enrichment)

# I think the sample is simply too small.

# Visualising findings
joined_pathways

dim(joined_pathways)
# 38  9

joined_pathways_tidy <- joined_pathways %>%
  separate(GeneRatio, into = c("numerator", "denominator")) %>%
  mutate(gene_ratio_nr = as.double(numerator) / as.double(denominator))

save(joined_pathways_tidy,file="joined_pathways_tidy.Rda")


pdf("pathway_enrichment_plot.pdf")
ggplot(data = joined_pathways_tidy) +
  geom_point(mapping = aes(x = gene_ratio_nr, y = Description, color = p.adjust, size = p.adjust)) +
  labs(
    title = paste(
      "Dotplot of pathway description against gene ratio"))
dev.off()

pdf("pathway_enrichment_plot_id.pdf")
ggplot(data = joined_pathways_tidy) +
  geom_point(mapping = aes(x = gene_ratio_nr, y = ID, color = p.adjust, size = p.adjust)) +
  labs(
    title = paste(
      "Dotplot of pathway ID against gene ratio"))
dev.off()

# Kegg analysis
atac_seq_control_treatment_significant
rna_seq_control_treatment_significant

save(atac_seq_control_treatment_significant,file="atac_seq_control_treatment_significant.Rda")
save(rna_seq_control_treatment_significant,file="rna_seq_control_treatment_significant.Rda")

TxDb.Hsapiens.UCSC.hg38.knownGene
files = c("atac_seq_control_treatment_significant", "rna_seq_control_treatment_significant")

peakAnnoList <- lapply(files, annotatePeak, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,tssRegion=c(-3000, 3000), verbose=FALSE)


# gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
# pathway2 <- enrichPathway(gene)


# combined kegg pathway analysis
# ?compareCluster
# 
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# files <- getSampleFiles()
# files
# 
# files2 <- list(atac = "C:/Users/Tamara/Documents/data_challenge/data (1)/data/atac_seq_control_treatment_significant.Rda", rna = "C:/Users/Tamara/Documents/data_challenge/data (1)/data/rna_seq_control_treatment_significant.Rda")
# files2
# file.choose()
# 
# 
# peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# peakAnnoList
# 
# peakAnnoList2 <- lapply(files2, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
# peakAnnoList2
# 
# genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) = sub("_", "/n", names(genes))
# summary(genes)
# 
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            pvalueCutoff  = 0.05,
#                            pAdjustMethod = "BH")
# 
# dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

# Peak annotation comparison with Kegg
anno_atac_peaks
anno_rna_peaks

anno_atac_peaks_df$geneId


atac_enrichment <- enrichKEGG(anno_atac_peaks_df$geneId, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)

pdf("Kegg_ATACseq_pathway_enrichment.pdf")
dotplot(atac_enrichment, showCategory = 15, title = "KEGG ATAC-seq Pathway Enrichment Analysis")
dev.off()

rna_enrichment <- enrichKEGG(anno_rna_peaks_df$geneId, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)

pdf("Kegg_RNA_seq_pathway_enrichment.pdf")
dotplot(rna_enrichment, showCategory = 15, title = "KEGG RNA-seq Pathway Enrichment Analysis")
dev.off()


# Joining the dataframes
atac_enrichment_df <- as.data.frame(atac_enrichment)
rna_enrichment_df <- as.data.frame(rna_enrichment)

# Identifying pathways in both programs
Kegg_joined_pathways <- semi_join(atac_enrichment_df, rna_enrichment_df, by = "ID")

dim(Kegg_joined_pathways)
# 84  9
# Kegg identifies more pathways

Kegg_joined_pathways_tidy <- Kegg_joined_pathways %>%
  arrange(desc(p.adjust)) %>%
  separate(GeneRatio, into = c("numerator", "denominator")) %>%
  mutate(gene_ratio_nr = as.double(numerator) / as.double(denominator)) 

#%>%
#  mutate(Description_fct = as.factor(Description))
# ?desc
# 
# order <- order(Kegg_joined_pathways_tidy$p.adjust, decreasing = FALSE)
# 
# order(Kegg_joined_pathways_tidy$p.adjust, )
# ordered_Kegg <- Kegg_joined_pathways_tidy[order, ]
# ordered_Kegg
# 
# desc(Kegg_joined_pathways_tidy$p.adjust)
# 
# Kegg_joined_pathways$Description

save(Kegg_joined_pathways_tidy,file="Kegg_joined_pathways_tidy.Rda")


pdf("Kegg_pathway_enrichment_plot.pdf")
ggplot(data = Kegg_joined_pathways_tidy) +
  geom_point(mapping = aes(x = gene_ratio_nr, y = Description, color = p.adjust, size = p.adjust)) +
  labs(
    title = paste(
      "Dotplot of Kegg pathway description against gene ratio"))
dev.off()

pdf("Kegg_pathway_enrichment_plot_id.pdf")
ggplot(data = Kegg_joined_pathways_tidy) +
  geom_point(mapping = aes(x = gene_ratio_nr, y = ID, color = p.adjust, size = p.adjust)) +
  labs(
    title = paste(
      "Dotplot of Kegg pathway ID against gene ratio"))
dev.off()


# Gene coverage
ATAC_cov <- covplot(myPeaks, xlab = "Chromosome Size (bp)",
                    ylab = "Chromosome ID",
                    title = "ATAC-seq Peaks over Chromosomes")
RNA_cov <- covplot(myPeaks_rna_seq, xlab = "Chromosome Size (bp)",
                   ylab = "Chromosome ID",
                   title = "RNA-seq Peaks over Chromosomes")


pdf("ATAC_peak_coverage_plot.pdf")
covplot(myPeaks, xlab = "Chromosome Size (bp)", ylab = "Chromosome ID", title = "ATAC-seq Peaks over Chromosomes")
dev.off()
# The function works, but the issue might lie in the nr of peaks or perhaps chromosomal assignment is incertain due to shorter fragment length.

pdf("RNA_peak_coverage_plot.pdf")
covplot(myPeaks_rna_seq, xlab = "Chromosome Size (bp)", ylab = "Chromosome ID", title = "RNA-seq Peaks over Chromosomes")
dev.off()

# Coverage plot on consensus peaks
joined_df_geneId
dim(joined_df_geneId)
# 446  20

joined_df_geneId_tidy <- joined_df_geneId %>%
  unite("peakid", start:end, sep = "-", remove = FALSE) 


GR_consensus_peaks <- GRanges(seqnames = joined_df_geneId_tidy$seqnames,
                              ranges = joined_df_geneId_tidy$peakid,)

?covplot

GR_consensus_peaks
covplot(GR_consensus_peaks)

pdf("Coverage_plot_RNA_and_ATAC_seq_consensus_peaks.pdf")
covplot(GR_consensus_peaks, xlab = "Chromosome Size (bp)",
        ylab = "Chromosome ID",
        title = "Consensus Peaks of RNA- and ATAC-seq over Chromosomes")
dev.off()

# Saving the joined IDs
Kegg_joined_IDs <- Kegg_joined_pathways_tidy$ID
Kegg_joined_IDs

# "hsa05132" "hsa04914" "hsa05017" "hsa05219" "hsa04670" "hsa04211" "hsa04370" "hsa04213" "hsa04668" "hsa05203"
# "hsa04910" "hsa05211" "hsa04662" "hsa04144" "hsa05224" "hsa05167" "hsa04931" "hsa04810" "hsa05160" "hsa04210"
# "hsa05170" "hsa05162" "hsa04218" "hsa04064" "hsa04625" "hsa05166" "hsa05417" "hsa05131" "hsa04917" "hsa05213"
# "hsa04919" "hsa05231" "hsa05163" "hsa00230" "hsa05218" "hsa04666" "hsa05215" "hsa04664" "hsa04611" "hsa04658"
# "hsa04140" "hsa05161" "hsa05202" "hsa04926" "hsa04722" "hsa04935" "hsa04912" "hsa04062" "hsa04930" "hsa04520"
# "hsa05100" "hsa05225" "hsa05235" "hsa05210" "hsa04068" "hsa05165" "hsa01522" "hsa05221" "hsa05223" "hsa04928"
# "hsa05212" "hsa04750" "hsa04380" "hsa04659" "hsa04725" "hsa04071" "hsa00562" "hsa04010" "hsa05135" "hsa05220"
# "hsa05205" "hsa04934" "hsa04072" "hsa04070" "hsa05214" "hsa04390" "hsa04360" "hsa04012" "hsa05142" "hsa04015"
# "hsa04510" "hsa04933" "hsa05222" "hsa04660"




getwd()

# Question 3: Identifying candidate TFs for regulation
BiocManager::install("enrichTF")
library(enrichTF)

# files
makebedtable(atac_seq_control_treatment_significant, "atac_seq_control_treatment_significant.html", basedirectory = getwd())

makebedtable(rna_seq_control_treatment_significant, "rna_seq_control_treatment_significant.html", basedirectory = getwd())
atacseq_peaks_bed <- read.delim(file = "atacseq_peaks.bed", header = FALSE, col.names = c("chr", "start", "end"))
atac_seq_control_treatment_significant
?read.delim

head(atacseq_peaks_bed)
# chr  start    end
# 1   1  10002  10507
# 2   1  20220  22634
# 3   1  28573  30038
# 4   1  37946  39588
# 5   1  90823  91498
# 6   1 107005 107256

system.file(package = "enrichTF", )

save(atacseq_peaks_bed,file="atacseq_peaks_bed.bed")
# The format might require Chr prefix for column 1.

foregroundBedPath <- system.file(package = "enrichTF", "vignettes","atacseq_peaks_bed.bed")
# The path may be wrong
foregroundBedPath

# data.file <- "C:\Users\Tamara\Documents\data_challenge\data (1)"
file.path("C:", "Users", "Tamara", "Documents", "data_challenge", "data_(1)", "data", "atacseq_peaks_bed.bed")
path_to_file <- file.path("C:", "Users", "Tamara", "Documents", "data_challenge", "data_(1)", "data", "atacseq_peaks_bed.bed")
path_to_file
getwd()

BiocManager::install("all_motif_rmdup")
library(all_motif_rmdup)

foregroundBedPath <- system.file(package = "enrichTF", system.file(path_to_file))
foregroundBedPath

gen <- genBackground(inputForegroundBed = foregroundBedPath)

# PECA_TF_enrich(inputForegroundBed = foregroundBedPath, genome = "hg38") 
# 
# PECA_TF_enrich()
# hg38 will have new functions

# ?PECA_TF_enrich
# PECA_TF_enrich(inputForegroundBed, genome, threads = 2, pipeName = "pipe", ...)
printMap()

setGenome("testgenome")

sampleAnnotFn <- file.path("tcells", "samples.tsv")

foregroundBedPath <- file.path("atac_bed","atacseq_peaks_bed.bed")

foregroundBedPath <- system.file(package = "enrichTF", "atac_bed","atacseq_peaks_bed.bed")


foregroundBedPath <- system.file(package = "enrichTF", "extdata","atacseq_peaks_bed")
gen <- genBackground(inputForegroundBed = foregroundBedPath)

?genBackground
# inputForegroundBed should be a character scalar, the directory of the foreground BED file


# Trying a different function
# BiocManager::install("TFutils")
# BiocManager::install("GenomicFiles")
# 
# 
# TFutils::tftColl
# 
# grep("NFK", names(TFutils::tftColl), value=TRUE)
# 
# #
# library(GenomicFiles)
# data(fimo16)
# fimo16

# Trying different functions
# install dependencies
devtools::install_github("demuellae/muLogR")
devtools::install_github("demuellae/muRtools")
devtools::install_github("demuellae/muReportR")

# install ChrAccR
devtools::install_github("GreenleafLab/ChrAccR", dependencies=TRUE)

# hg38 annotation package
install.packages("https://muellerf.s3.amazonaws.com/data/ChrAccR/data/annotation/ChrAccRAnnotationHg38_0.0.1.tar.gz")
library(ChrAccR)

devtools::install_github("GreenleafLab/ChrAccRex")

# visualising practice dataset
library(ggplot2)
# use a grid-less theme
theme_set(muRtools::theme_nogrid())

##############################
# Retrying old code
#sampleAnnotFn <- file.path("tcells", "samples.tsv")
library(enrichTF)

foregroundBedPath <- file.path("atac_bed","atacseq_peaks_bed.Rda")
foregroundBedPath
# "atac_bed/atacseq_peaks_bed.bed"

data.file <- system.file(foregroundBedPath, package =  "enrichTF", lib.loc = NULL, mustWork = FALSE)
# Could I use other packages?

foregroundBedPath <- system.file(package = "enrichTF", "atac_bed","atacseq_peaks_bed.bed")

gen <- genBackground(inputForegroundBed = foregroundBedPath)

# Install 
MotifDb::associateTranscriptionFactors()
?associateTranscriptionFactors

MotifDb::geneToMotif()

# Peaks
motif_peaks <- read.table("atacseq_peaks.bed", header=FALSE)

atac_seq_control_treatment_significant
atac_seq_control_treatment_significant_df <- as.data.frame(atac_seq_control_treatment_significant)



atac_seq_control_treatment_significant_df <- atac_seq_control_treatment_significant_df[order(atac_seq_control_treatment_significant_df$padj)]

atac_seq_control_treatment_significant <- atac_seq_control_treatment_significant[order(atac_seq_control_treatment_significant$padj)]
atac_peaks_resized <-  resize(atac_seq_control_treatment_significant, width = 50, fix='center')

install.packages("rGADEM")
BiocManager::install("rGADEM")
library(rGADEM)



motif_peaks <- makeGRangesFromDataFrame(motif_peaks, keep.extra.columns = TRUE)

