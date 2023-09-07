# Txi subsetting

# load libaraies
library(tidyverse, warn.conflicts = FALSE)
library(tximport)
library(DESeq2)


# load data
load(file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/sample-gene-lists-forTXI.RData", verbose =TRUE)

load("/scratch/group/hu-lab/pacocean-metaT/tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)


# TXI subset function
## Subset txi directly
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

#
#
## Subset eukaryotes only, and keep all samples.
#
#
txi_euk_annot <- subsetTxi(txi, all_samples, euks_only)

## Subset
tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(all_samples$sample))
## Set names
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_euk_annot$counts)


ds_tpm_light <- DESeqDataSetFromTximport(txi_euk_annot,
                                         colData = tmp_sample_merged,
                                         design = ~0 + LIGHT)

ds_tpm_pac <- DESeqDataSetFromTximport(txi_euk_annot,
                                       colData = tmp_sample_merged,
                                       design = ~0 + PACIFIC_REGION)

ds_tpm_lightpac <- DESeqDataSetFromTximport(txi_euk_annot,
                                            colData = tmp_sample_merged,
                                            design = ~0 + PACIFIC_REGION + LIGHT)

save(ds_tpm_light, ds_tpm_pac, ds_tpm_lightpac, file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/wholecomm-light-pacregion-deseq.RData")


