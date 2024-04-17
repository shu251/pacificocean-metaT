library(DESeq2)
library(tidyverse)
library(tximport)

# Import R objects we need
# load("/scratch/group/hu-lab/pacocean-metaT/tximport_metaT-ALOHA-CA.Rdata", verbose = TRUE)
# taxfxn <- read.csv("/scratch/group/hu-lab/pacocean-metaT/TaxonomicAndFunctionalAnnotations.csv")

# Isolate what is needed from tax and fxn annotations for most analyses.
# taxfxn_mini <- taxfxn %>%
  # select(SequenceID, Domain:Species, PFAMs, KEGG_ko)

# save(sample_merged, taxfxn_mini, file = "/scratch/group/hu-lab/pacocean-metaT/metadata_files_grace.RData")

load(file = "/home/skhu/pacificocean-metaT/input-data/annot_files.RData")

load(file = "/scratch/group/hu-lab/pacocean-metaT/normed_dfs_vst.RData")
# vst_all, vst_subset, df_ctr_norm_vst, df_ctr_norm_subset_vst

#load("/scratch/group/hu-lab/pacocean-metaT/normed_dfs_rld.RData")
# rld_all, rld_subset, df_ctr_norm_rld, df_ctr_norm_subset_rld

# function to convert into a usable transformed dataset
convert_format_scaled <- function(df){
  all_samples_col_header <- colnames(df)
  df %>% 
    # sample_n(1000) %>% #uncomment for testing
    rownames_to_column(var = "SequenceID") %>% 
    pivot_longer(cols = all_of(all_samples_col_header), names_to = "Sample_rep", values_to = "VALUE") %>% 
    left_join(sample_merged) %>% 
    unite(SAMPLE_NOREP, REGION, DEPTH_CATEGORY, sep = " ", remove = FALSE) %>% 
    group_by(SequenceID, SAMPLE_NOREP, REGION, DEPTH_CATEGORY, DEPTH, PACIFIC_REGION, LIGHT) %>% 
    summarize(MEAN_VALUE = mean(VALUE)) %>% 
    left_join(taxfxn_mini)
}

cat("\n\n Start fxn to compile output \n\n")

# Input dataframe
# rlog_format_subset <- convert_format_scaled(df_ctr_norm_subset_rld)
# rlog_format_all <- convert_format_scaled(df_ctr_norm_rld)

# save(rlog_format_subset, rlog_format_all, file = "/scratch/group/hu-lab/pacocean-metaT/rlog_formatted_forlocaluse.RData")
# cat("\n\n Done with rlog \n\n")
## ADD in toy data set

vst_format_subset <- convert_format_scaled(df_ctr_norm_subset_vst)
vst_format_all <- convert_format_scaled(df_ctr_norm_vst)

save(vst_format_subset, vst_format_all, file = "/scratch/group/hu-lab/pacocean-metaT/vst_formatted_forlocaluse.RData")

cat("\n\n Complete. \n\n")




