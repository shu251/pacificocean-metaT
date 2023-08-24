# Setting up for repeat DESeq - functions for running these comparisons.

library(DESeq2)
library(tidyverse)
library(tximport)


# Load data
load("tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)
# taxfxn <- read.csv("TaxonomicAndFunctionalAnnotations.csv")
loca("sample-gene-lists_TXISUBSET.RData", verbose = TRUE)

# Subset txi directly
subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
  {
  genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
  txi$abundance <- txi$abundance[genes, samples$sample]
  txi$counts <- txi$counts[genes, samples$sample]
  txi$length <- txi$length[genes, samples$sample]
  return(txi)
  }

# Example usage 
# tmp <- sample(taxfxn$SequenceID,10,replace = FALSE)
# pola <- data.frame(sample = c("PortofLA_1", "PortofLA_2"))
# txi_pola <- subsetTxi(txi, pola, tmp)

# Subset txi
txi_npsg <- subsetTxi(txi, npsg_only, genes_tax_fxn_all)
txi_ca <- subsetTxi(txi, ca_only, genes_tax_fxn_all)

save(txi_npsg, txi_ca, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/station_txi.RData")

ds_tpm_light <- DESeqDataSetFromTximport(txi,
                                         colData = sample_merged,
                                         design = ~0 + LIGHT)
save(ds_tpm_light, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/ds_tpm_light.RData")

ds_tpm_npsg_light <- DESeqDataSetFromTximport(txi_npsg,
                                              colData = sample_merged,
                                              design = ~0 + LIGHT)

ds_tpm_npsg_month <- DESeqDataSetFromTximport(txi_npsg,
                                              colData = sample_merged,
                                              design = ~0 + MONTH)

save(ds_tpm_npsg_light, ds_tpm_npsg_month, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/npsg_only.RData")

ds_tpm_ca_light <- DESeqDataSetFromTximport(txi_ca,
                                            colData = sample_merged,
                                            design = ~0 + LIGHT)
save(ds_tpm_ca_light, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/ds_tpm_ca_light.RData")

txi_euph <- subsetTxi(txi, euphotic, genes_tax_fxn_all)
txi_subeuph <- subsetTxi(txi, subeuphotic, genes_tax_fxn_all)
save(txi_euph, txi_subeuph, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/light_txi.RData")


ds_tpm_euphotic <- DESeqDataSetFromTximport(txi_euph,
                                            colData = sample_merged,
                                            design = ~0 + PACIFIC_REGION)

ds_tpm_subeuphotic <- DESeqDataSetFromTximport(txi_subeuph,
                                               colData = sample_merged,
                                               design = ~0 + PACIFIC_REGION)
save(ds_tpm_subeuphotic, ds_tpm_euphotic, file = "/vortexfs1/scratch/sarahhu/txi-objs-metaT/light_deseqs.RData")


# counts_scaled <- makeCountsFromAbundance(
#   as.matrix(txi$counts),
#   as.matrix(txi$abundance),
#   as.matrix(txi$length),
#   countsFromAbundance = "scaledTPM"
# )
# 
# counts_df <- as.data.frame(counts_scaled)
# 
# 
# ds_tpm_samplename <- DESeqDataSetFromTximport(txi,
#                                               colData = sample_merged,
#                                               design = ~0 + SAMPLENAME)














