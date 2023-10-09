# load libaraies
library(tidyverse, warn.conflicts = FALSE)
library(tximport)
library(DESeq2)

# Load data
load(file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/sample-gene-lists-forTXI.RData", verbose =TRUE)
load("/scratch/group/hu-lab/pacocean-metaT/tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)

# Function to subset txi by individual taxa
deseq_region_bytax <- function(sample_set, gene_set){
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
  # Import as DESeq object - use PACIFIC_REGION in the design
  # DESeq
  ds_tpm_output <- DESeqDataSetFromTximport(txi_output,
                                            colData = tmp_sample_merged,
                                            design = ~0 + PACIFIC_REGION)
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
  ## Positive log fold change == up regulared in CA, compared to NPSG
  ds_tpm_output_filtered$PACIFIC_REGION <- factor(ds_tpm_output_filtered$PACIFIC_REGION, levels = c("NPSG", "CA"))
  de_tax_output <- DESeq2::DESeq(ds_tpm_output_filtered)
  resultsNames(de_tax_output)
  summary(de_tax_output)
  cat("\n\nCompleted.\n\n")
  return(de_tax_output)
}

## Apply to each taxa

# From the euphotic samples, what are DE genes among haptophytes between CA and NPSG?
cat("\n\nDiatoms start.\n\n")
de_diatom <- deseq_region_bytax(euphotic, diatom)


cat("\n\nHaptophytes start.\n\n")
# From the euphotic samples, what are DE genes among haptophytes between CA and NPSG?
de_hapto <- deseq_region_bytax(euphotic, hapto)


cat("\n\nDinoflagellates start.\n\n")
# From the euphotic samples, what are DE genes among Dinoflagellates between CA and NPSG?
de_dinos <- deseq_region_bytax(euphotic, dinos)

cat("\n\nChlorophytes start.\n\n")
# From the euphotic samples, what are DE genes among chlorophytes between CA and NPSG?
de_chloro <- deseq_region_bytax(euphotic, chloro)

save(de_diatom, de_hapto, de_dino, de_chloro, file = "/scratch/group/hu-lab/pacocean-metaT/Robjs/euphotic_by_taxa.RData")