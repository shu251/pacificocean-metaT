# Set up DESeq comparison parameters

## Library loading
library(tidyverse, warn.conflicts = FALSE)

## Import needed Robjects
metadata <- read.csv("sample-list-revised.txt")
taxfxn <- read.csv("TaxonomicAndFunctionalAnnotations.csv")
key_geneid <- read.csv("input-data/keygene_id.csv")

###############
### SAMPLES ###
###############

all_samples <- metadata %>% 
  select(sample = SAMPLE)

npsg_only <- metadata %>% 
  filter(PACIFIC_REGION == "NPSG") %>% 
  select(sample = SAMPLE)
# Compare euphotic and subeuphotic
# Use this to compare across months too

ca_only <- metadata %>% 
  filter(PACIFIC_REGION != "NPSG") %>% 
  select(sample = SAMPLE)
# Compare euphotic and subeuphotic

euphotic <- metadata %>% 
  filter(LIGHT == "Euphotic")%>% 
  select(sample = SAMPLE)
# Compare NPSG to CA euphotic zone

subeuphotic <- metadata %>% 
  filter(LIGHT != "Euphotic")%>% 
  select(sample = SAMPLE)
# Compare NPSG to CA subeuphotic zone

############
## Genes ###
############

# Get list of annotations where a GO, PFAM, and KEGG ID were assigned.
genes_fxn_all <- as.character(
  filter(taxfxn, GOs != "-" & PFAMs != "-" & KEGG_ko != "-") %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]]) #this line outputs the selected vector from the pipe

# Removing unclassified and unannotated sequences
genes_tax_fxn_all <- as.character(
  filter(taxfxn, GOs != "-" & PFAMs != "-" & KEGG_ko != "-"
         & Domain == "Eukaryota" & Supergroup != "Unclassified") %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])


kegg_keep <- as.character(key_geneid %>% 
                            filter(Category01 != "" & !(is.na(Category01))) %>% 
                            select(KEGG) %>% 
                            .[["KEGG"]])
# length(kegg_keep)

genes_kegg_curated <- as.character(taxfxn %>%
                                     mutate(KEGG_mod = str_remove_all(KEGG_ko, "ko:")) %>% 
                                     filter(KEGG_mod %in% kegg_keep) %>% 
                                     select(SequenceID) %>% 
                                     .[["SequenceID"]])

############
### TAXA ###
############
euks_only <- as.character(
  filter(taxfxn, Domain == "Eukaryota" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])

diatom <- as.character(
  filter(taxfxn, Class == "Bacillariophyta" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]]) 

#Dinoflagellata Phylum
dinos <- as.character(
  filter(taxfxn, Phylum == "Dinoflagellata" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]]) 

#Haptophyta in phylum
hapto <- as.character(
  filter(taxfxn, Phylum == "Haptophyta" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])

# Chlorophyta
chloro <- as.character(
  filter(taxfxn, Phylum == "Chlorophyta" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])

#Ciliophora
ciliate <- as.character(
  filter(taxfxn, Phylum == "Ciliophora" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])

# Rhizaria at supergroup
rhizaria <- as.character(
  filter(taxfxn, Supergroup == "Rhizaria" & (PFAMs != "-" | KEGG_ko != "-")) %>% 
    select(SequenceID) %>% 
    .[["SequenceID"]])

############
# Save all #
############

save(
  # These are data frames
  all_samples, npsg_only, ca_only, euphotic, subeuphotic,
  # these are lists (character lists)
  genes_fxn_all, genes_tax_fxn_all, key_geneid, genes_kegg_curated, 
  euks_only, diatom, dinos, hapto, chloro, ciliate, rhizaria,
  file = "sample-gene-lists-forTXI.RData")