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
# By region
#
#
# NPSG
#
txi_npsg <- subsetTxi(txi, npsg_only, euks_only)


## Subset
tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(npsg_only$sample))
## Set names
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_npsg$counts)

# Compare Euphotic vs. sub-euphotic samples in the NPSG
ds_tpm_npsg_light <- DESeqDataSetFromTximport(txi_npsg,
                                              colData = tmp_sample_merged,
                                              design = ~0 + LIGHT)

# Compare July vs. March in the NPSG
ds_tpm_npsg_month <- DESeqDataSetFromTximport(txi_npsg,
                                              colData = tmp_sample_merged,
                                              design = ~0 + MONTH)

save(ds_tpm_npsg_light, ds_tpm_npsg_month, file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/npsg-deseq.RData")


# 
#
#
# CA
#
#
txi_ca <- subsetTxi(txi, ca_only, euks_only)

tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(ca_only$sample))
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_ca$counts)

# Compare euphotic vs. subeuphotic in coastal California
## Includes Port of LA and Catalina
ds_tpm_ca_light <- DESeqDataSetFromTximport(txi_ca,
                                            colData = tmp_sample_merged,
                                            design = ~0 + LIGHT)

save(ds_tpm_ca_light, file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/ca-deseq.RData")
