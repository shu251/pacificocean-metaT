setwd("/scratch/user/skhu/pacocean")
library(DESeq2)
library(tidyverse)
library(tximport)

load(file = "sample-gene-lists-forTXI.RData", verbose =TRUE)
load(file = "tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)

# Isolating dinoflagellates, ciliates, rhizaria, and haptophytes

# Function to subset txi by individual taxa
deseq_depth_bytax <- function(sample_set, gene_set){
  # First incorporate the txi subset fxn
  subsetTxi <- function(txi, samples, include_genes=rownames(txi$counts))
  {
    genes <- rownames(txi$counts)[rownames(txi$counts) %in% include_genes]
    txi$abundance <- txi$abundance[genes, samples$sample]
    txi$counts <- txi$counts[genes, samples$sample]
    txi$length <- txi$length[genes, samples$sample]
    return(txi)
  }
  # Run subset and sample merge re-set
  txi_output <- subsetTxi(txi, sample_set, gene_set)
  # This script will keep replacing tmp_sample_merged
  tmp_sample_merged <- sample_merged %>% 
    filter(Sample_rep %in% as.character(sample_set$sample))
  rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
  rownames(tmp_sample_merged) <- colnames(txi_output$counts)
  #
  # Import as DESeq object - use LIGHT in the design
  # DESeq
  ds_tpm_output <- DESeqDataSetFromTximport(txi_output,
                                            colData = tmp_sample_merged,
                                            design = ~0 + LIGHT)
  # return(ds_tpm_output)
  # Further process DESeq
  groupsize <- 2 # Transcript to consider, must be in at least 3 samples
  keep <- rowSums(counts(ds_tpm_output) >= 10) >= groupsize # And have >= to 10 counts
  ds_tpm_output_filtered <- ds_tpm_output[keep,]
  ###
  # Filtering stats:
  cat("\nStarted with ", dim(ds_tpm_output)[1], "observations. Filtering by 2 samples and 10 counts resulted in,", dim(ds_tpm_output_filtered)[1], ", which is", (100*(dim(ds_tpm_output_filtered)[1]/dim(ds_tpm_output)[1])), "% of the data.\n\n")
  ###
  #
  ## Positive log fold change == up regularted in Euphotic, compared to Sub-euphotic
  ds_tpm_output_filtered$LIGHT <- factor(ds_tpm_output_filtered$LIGHT, levels = c("Sub-euphotic", "Euphotic"))
  de_tax_output <- DESeq2::DESeq(ds_tpm_output_filtered)
  resultsNames(de_tax_output)
  summary(de_tax_output)
  cat("\n\nCompleted.\n\n")
  return(de_tax_output)
}

## Apply to each taxa

# From the all samples, what are the differences in euphotic vs sub-eiphotic?
cat("\n\nStart with all samples.\n\n")
cat("\n\nCiliates.\n\n")
de_all_ciliate <- deseq_depth_bytax(all_samples, ciliate)
cat("\n\nDinos.\n\n")
de_all_dino <- deseq_depth_bytax(all_samples, dinos)
cat("\n\nHaptophytes.\n\n")
de_all_hapto <- deseq_depth_bytax(all_samples, hapto)
cat("\n\nRhizaria.\n\n")
de_all_rhiz <- deseq_depth_bytax(all_samples, rhizaria)

save(de_all_ciliate, de_all_dino, de_all_hapto, de_all_rhiz, file = "/scratch/user/skhu/pacocean/ALL_LIGHT_by_taxa.RData")

# From the NPSG samples, what are the differences in euphotic vs sub-eiphotic?
cat("\n\nMovig to NPSG samples.\n\n")
cat("\n\nCiliates.\n\n")
de_npsg_ciliate <- deseq_depth_bytax(npsg_only, ciliate)
cat("\n\nDinos.\n\n")
de_npsg_dino <- deseq_depth_bytax(npsg_only, dinos)
cat("\n\nHaptophytes.\n\n")
de_npsg_hapto <- deseq_depth_bytax(npsg_only, hapto)
cat("\n\nRhizaria.\n\n")
de_npsg_rhiz <- deseq_depth_bytax(npsg_only, rhizaria)

save(de_npsg_ciliate, de_npsg_dino, de_npsg_hapto, de_npsg_rhiz, file = "/scratch/user/skhu/pacocean/NPSG_LIGHT_by_taxa.RData")


# From the NPSG samples, what are the differences in euphotic vs sub-eiphotic?
cat("\n\nMovig to CA samples.\n\n")
cat("\n\nCiliates.\n\n")
de_ca_ciliate <- deseq_depth_bytax(ca_only, ciliate)
cat("\n\nDinos.\n\n")
de_ca_dino <- deseq_depth_bytax(ca_only, dinos)
cat("\n\nHaptophytes.\n\n")
de_ca_hapto <- deseq_depth_bytax(ca_only, hapto)
cat("\n\nRhizaria.\n\n")
de_ca_rhiz <- deseq_depth_bytax(ca_only, rhizaria)

save(de_ca_ciliate, de_ca_dino, de_ca_hapto, de_ca_rhiz, file = "/scratch/user/skhu/pacocean/CA_LIGHT_by_taxa.RData")