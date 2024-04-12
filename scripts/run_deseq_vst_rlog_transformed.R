setwd("/home/skhu/pacificocean-metaT/scripts")
library(DESeq2)
library(tidyverse)
library(tximport)

# Import R objects we need
load("/scratch/group/hu-lab/pacocean-metaT/tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)

taxfxn <- read.csv("/scratch/group/hu-lab/pacocean-metaT/TaxonomicAndFunctionalAnnotations.csv")

# Isolate what is needed from tax and fxn annotations for most analyses.
taxfxn_mini <- taxfxn %>%
  select(SequenceID, Domain:Species, PFAMs, KEGG_ko)

# Make a DESeq object
ds_tpm_samplename <- DESeqDataSetFromTximport(txi,
                                              colData = sample_merged,
                                              design = ~0 + SAMPLENAME)

## Create a subsampled dataset - remove low count transcripts
keep <- rowSums(counts(ds_tpm_samplename)) >= 10
dds_subset <- ds_tpm_samplename[keep,]



# Run DESeq
## Estimates dispersion\

dds_dis_all <- DESeq(ds_tpm_samplename)
# res_all <- results(dds_dis_all)

dds_dis_subset <- DESeq(dds_subset)
# res_subset <- results(dds_dis_subset)

cat("\n\n Done with DESeq steps \n\n")

# Estiamte scaled datasets by VST
vst_all <- vst(dds_dis_all, blind = FALSE) 
vst_subset <- vst(dds_dis_subset, blind = FALSE)

# Get data frames 
df_ctr_norm_vst <- as.data.frame(assay(vst_all))
df_ctr_norm_subset_vst <- as.data.frame(assay(vst_subset))

save(vst_all, vst_subset, df_ctr_norm_vst, df_ctr_norm_subset_vst,
     file = "/scratch/group/hu-lab/pacocean-metaT/normed_dfs_vst.RData")

cat("\n\n Saved VST output \n\n")

# Repeat for rlog transformed data
rld_all <- rlog(dds_dis_all, blind = FALSE)
rld_subset <- rlog(dds_dis_subset, blind = FALSE)

df_ctr_norm_rld <- as.data.frame(assay(rld_all))
df_ctr_norm_subset_rld <- as.data.frame(assay(rld_subset))

save(rld_all, rld_subset, df_ctr_norm_rld, df_ctr_norm_subset_rld, file = "/scratch/group/hu-lab/pacocean-metaT/normed_dfs_rld.RData")
cat("\n\n Saved Rlog output \n\n")
