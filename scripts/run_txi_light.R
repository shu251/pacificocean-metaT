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
#
# light availability only
#
#

txi_euph <- subsetTxi(txi, euphotic, euks_only)

tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(euphotic$sample))
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_euph$counts)

# Compare CA versus NPSG within euphotic samples
ds_tpm_euphotic <- DESeqDataSetFromTximport(txi_euph,
                                            colData = tmp_sample_merged,
                                            design = ~0 + PACIFIC_REGION)



txi_subeuph <- subsetTxi(txi, subeuphotic, euks_only)

tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(subeuphotic$sample))
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_subeuph$counts)

# Compare CA versus NPSG at subeuphotic depths
ds_tpm_subeuphotic <- DESeqDataSetFromTximport(txi_subeuph,
                                               colData = tmp_sample_merged,
                                               design = ~0 + PACIFIC_REGION)

save(ds_tpm_subeuphotic, ds_tpm_euphotic, file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/euphotic-subeuphotic-deseq.RData")
