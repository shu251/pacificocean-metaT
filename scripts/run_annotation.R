library(tidyverse)

# Import R objects we need

# load(file = "/home/skhu/pacificocean-metaT/input-data/annot_files.RData")
# 
# load("/scratch/group/hu-lab/pacocean-metaT/TPM-count-tables.RData")
# 
# mean_counts_annotated <- mean_counts_df %>% 
#   rownames_to_column(var = "SequenceID") %>%
#   left_join(taxfxn_mini) %>%
#   pivot_longer(cols = starts_with("mean."), names_to = "sample_tmp", values_to = "TPM") %>% 
#   filter(TPM > 0) %>% 
#   mutate(Sample = str_remove(sample_tmp, "mean.")) %>% 
#   left_join(sample_merged) %>% 
#   select(-sample_tmp, -SAMPLENAME) %>% 
#   unite(SAMPLENAME, PACIFIC_REGION, REGION, DEPTH_CATEGORY, sep = " ", remove = FALSE)
# 
# save(mean_counts_annotated, file = "/scratch/group/hu-lab/pacocean-metaT/mean_counts_annotated.RData")
# cat("\n\n Complete. \n\n")
load("/scratch/group/hu-lab/pacocean-metaT/mean_counts_annotated.RData", verbose = TRUE)

# Import kegg
kegg <- read.csv("/home/skhu/KEGG_DB/combined_kegg.csv")
curated_kegg <- read.csv("/home/skhu/KEGG_DB/reformat-kegg-pfam-skh.csv")

key_geneid <- curated_kegg %>% 
  select(-X) %>% 
  right_join(kegg %>% select(KEGG = KO_number, everything(), -X)) %>% 
  distinct() %>% 
  select(starts_with("KeggOrthology_"), Category01, Category02, FullName, GeneID, Gene_identification, KEGG, PFAM, Descriptions, REF = REFs)

kegg_ortho_based <- key_geneid %>% 
  select(KeggOrthology_B, KEGG, GeneID, Gene_identification) %>% 
  filter(!is.na(KeggOrthology_B)) %>% 
  distinct()

kegg_curated <- key_geneid %>% 
  select(Category01, Category02, KEGG, GeneID, Gene_identification) %>% 
  filter(!is.na(Category01)) %>% 
  distinct()
# head(kegg_curated)

cat("\n\n A look at mean_counts_annotated\n\n")
glimpse(mean_counts_annotated)

#
cat("\n\n Total unique transcripts:\n")
length(unique(mean_counts_annotated$SequenceID))

cat("\n\n Sum TPM of everything:\n")
sum(mean_counts_annotated$TPM)

glimpse(key_geneid)
# 
cat("\n\n Isolate curated TPMs only\n\n")

mean_counts_curated_only <- curated_kegg %>% 
  left_join(mean_counts_annotated %>% 
              mutate(KEGG = str_remove(KEGG_ko, "ko:")))

save(mean_counts_curated_only, file = "/scratch/group/hu-lab/pacocean-metaT/mean_counts_curated_only.RData")
  
cat("\n\n A look at mean_counts_curated_only\n\n")
glimpse(mean_counts_curated_only)

#
cat("\n\n Total unique transcripts:\n")
length(unique(mean_counts_curated_only$SequenceID))

cat("\n\n Sum TPM of everything:\n")
sum(mean_counts_curated_only$TPM)
  
cat("\n\n Complete.\n")
  
  
  
