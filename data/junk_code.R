# JUnk




# rna_seq_control_treatment_significant <- rna_seq_control_treatment[(!is.na(rna_seq_control_treatment$padj) & 
#                                                                       rna_seq_control_treatment$padj < 0.05) & rna_seq_control_treatment %over% toOverLap_genes, ]
# 
# rna_seq_control_treatment_significant
# The output is the same, whether you are using genes or promoters for the annotation.



# and save the deseq results.
save(rna_seq_control_treatment_significant, file="significant_deseq2_rna-seq_peaks.Rda")

rna_seq_control_treatment_significant

significant_atac_peaks_df <- as.data.frame(atac_seq_control_treatment_significant)
dim(significant_atac_peaks_df)
# 17978    11

significant_rna_peaks_df <- as.data.frame(rna_seq_control_treatment_significant)
dim(significant_rna_peaks_df)
# 6363   11

significant_atac_peaks_df
significant_rna_peaks_df

# joining the start and end coordinates
significant_atac_peaks_df <- significant_atac_peaks_df %>%
  unite("peakid", start:end, sep = "-", remove = FALSE)

significant_atac_peaks_df


significant_rna_peaks_df <- significant_rna_peaks_df %>%
  unite("peakid", start:end, sep = "-", remove = FALSE)

significant_rna_peaks_df

significant_rna_peaks_df$peakid
significant_atac_peaks_df$peakid

summary(significant_rna_peaks_df$peakid)
# Length     Class      Mode 
# 6466 character character

summary(significant_atac_peaks_df$peakid)
# Length     Class      Mode 
# 29524 character character 

which(significant_rna_peaks_df$peakid %in% significant_atac_peaks_df$peakid)
# 0
# I need to subset those RNA seq peaks that are differentially open in the ATAC-seq data (among the treated group), 


significant_rna_peaks_df
significant_atac_peaks_df

rna_seq_control_treatment_significant
atac_seq_control_treatment_significant

significant_atac_peaks_df$start

consistent_peaks_atac_seq_subset <- subset(discretised_atac_seq_peaks_df, (discretised_atac_seq_peaks_df$s84 > 0 & atacseq_peak_counts$s85 > 0 &atacseq_peak_counts$s86 > 0) | (atacseq_peak_counts$s93 > 0 & atacseq_peak_counts$s94 > 0 &atacseq_peak_counts$s95 > 0))
dim(consistent_peaks_atac_seq_subset)
# 173284      6

consistent_peak_counts_atac_seq_subset <- subset(atacseq_peak_counts, (atacseq_peak_counts$s84 > 0 & atacseq_peak_counts$s85 > 0 &atacseq_peak_counts$s86 > 0) | (atacseq_peak_counts$s93 > 0 & atacseq_peak_counts$s94 > 0 &atacseq_peak_counts$s95 > 0))
dim(consistent_peak_counts_atac_seq_subset)
# 173284      8

# Subset those significant atac-seq peaks that overlap one or more of the significant RNA-seq peaks 
subset(significant_atac_peaks_df, (significant_atac_peaks_df) )

##############################################################
significant_rna_peaks_df
summary(consistent_peak_counts_atac_seq_subset$chr)


myPeaks <- GRanges(
  seqnames = consistent_peak_counts_atac_seq_subset$chr,
  ranges = consistent_peak_counts_atac_seq_subset$`start-end`,
  control_1 = consistent_peak_counts_atac_seq_subset$s84,
  control_2 = consistent_peak_counts_atac_seq_subset$s85,
  control_3 = consistent_peak_counts_atac_seq_subset$s86,
  treated_1 = consistent_peak_counts_atac_seq_subset$s93,
  treated_2 = consistent_peak_counts_atac_seq_subset$s94,
  treated_3 = consistent_peak_counts_atac_seq_subset$s95,)

myPeaks

# Subsetting those regions of open chromatin within promoters
# and creating a table to be able to view the results in IGV
# using the makebedtable function from the trachtables package
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# 
# toOverLap <- promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, 2000, 2000)
# toOverLap
# 
# 
toOverLap <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "gene_id", filter=NULL, single.strand.genes.only=TRUE)
toOverLap

atac_seq_control_treatment



# 
# toOverLap <- cds(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "cds_id", filter=NULL, use.names=FALSE)
# toOverLap
# ?genes
# 
# TxDb.Hsapiens.UCSC.hg38.knownGene

toOverLap
# subsetting regions of significantly enriched open chromatin
# atac_seq_control_treatment_significant_overlap <- atac_seq_control_treatment[(!is.na(atac_seq_control_treatment$padj) & 
#                                                                           atac_seq_control_treatment$padj < 0.05) & atac_seq_control_treatment %over% toOverLap, ]
# atac_seq_control_treatment_significant_overlap
# Warning message:
#   In .Seqinfo.mergexy(x, y) :
#   The 2 combined objects have no sequence levels in common. (Use
#                                                              suppressWarnings() to suppress this warning.)
# This command cannot be used for the atac or the RNA seq data


# makebedtable(atac_seq_control_treatment_significant, "atac_seq_control_treatment_significant.html")




# sig_atac_df <- as.data.frame(atac_seq_control_treatment_significant)
# view(sig_atac_df)
# 
#  
# sig_atac_df <- sig_atac_df %>%
#   unite("ranges", start:end, sep = ",") 
# sig_atac_df

# 
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
# Error in asMethod(object) : 
# The character vector to convert to an IRanges object must contain strings of the form
# "start-end" or "start..end", with end >= start - 1, or just "pos". For example: "2501-2900",
# "2501..2900", or "740".

# sig_rna_df <- as.data.frame(rna_seq_control_treatment_significant)
# view(sig_rna_df)


# I cannot map the atac-seq reads to human genome build 38. The warning message states that human genome build 38 is not accessible 
# for my version of R. I will try to debug this later.

# # Transcription start sites
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), fix = "start", 1)
# TSSs
# 
# 