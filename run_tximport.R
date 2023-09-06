# Last updated Aug 2023 - SKH
## Inputs: salmon quant file locations
# SRR ids that correspond to salmon files
# TaxonomicAndFunctionalAnnotations.csv and TaxonomicAndFunctionalAnnotations_names.csv files from eukrhythmic
# Metadata information: complete-sample-list.txt
# Memory to run this Rscript

# Load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(tximport)
library(readr)

num_threads <- getDTthreads()

# Get all salmon output files
files <- Sys.glob("/proj/omics/huber/shu/metaT-SPOT-ALOHA/2023-akrinos-salmonrerun/Hu2023-SalmonRemap/salmon_sample/*quant/quant.sf")


# Create empty dataframe
sample_table <- data.frame("SRR"=c("SRR6048900","SRR6048899",
                                "SRR6048898","SRR6048897",
                                "SRR6048896","SRR6048895",
                                "SRR6048894","SRR6048893",
                                "SRR6048892","SRR6048891",
                                "SRR6048902","SRR6048901",
                                "SRR11178183","SRR11178182",
                                "SRR11178173","SRR11178172",
                                "SRR11178171","SRR11178170",
                                "SRR11178169","SRR11178168",
                                "SRR11178167","SRR11178166",
                                "SRR11178181","SRR11178180",
                                "SRR5799332","SRR5799333",
                                "SRR5799340","SRR5799341",
                                "SRR5799343","SRR5799344",
                                "SRR5799342","SRR5799338",
                                "SRR5799339","SRR5799336",
                                "SRR5799337","SRR5799334",
                                "SRR5799335","SRR11178179",
                                "SRR11178178","SRR11178177",
                                "SRR11178176","SRR11178175",
                                "SRR11178174"),
                        "SAMN"=c("SAMN07647713","SAMN07647714",
                                 "SAMN07647715","SAMN07647716",
                                 "SAMN07647717","SAMN07647718",
                                 "SAMN07647719","SAMN07647720",
                                 "SAMN07647721","SAMN07647722",
                                 "SAMN07647723","SAMN07647724",
                                 "SAMN14206057","SAMN14206058",
                                 "SAMN14206059","SAMN14206060",
                                 "SAMN14206061","SAMN14206062",
                                 "SAMN14206063","SAMN14206064",
                                 "SAMN14206065","SAMN14206066",
                                 "SAMN14206067","SAMN14206068",
                                 "SAMN07269832","SAMN07269833",
                                 "SAMN07269834","SAMN07269835",
                                 "SAMN07269836","SAMN07269837",
                                 "SAMN07269838","SAMN07269826",
                                 "SAMN07269827","SAMN07269828",
                                 "SAMN07269829","SAMN07269830",
                                 "SAMN07269831","SAMN14206069",
                                 "SAMN14206070","SAMN14206071",
                                 "SAMN14206072","SAMN14206073",
                                 "SAMN14206074"),
                        "Sample"=c("July_1000m","July_150m","July_5m","July_5m",
                                   "July_DCM","July_DCM","March_1000m","March_150m","March_5m","March_5m",
                                   "March_DCM","March_DCM","Catalina","Catalina",
                                   "Catalina","Catalina","Catalina","Catalina",
                                   "PortofLA",
                                   "PortofLA",
                                   "PortofLA",
                                   "PortofLA",
                                   "PortofLA",
                                   "PortofLA","SPOT_150m","SPOT_150m","SPOT_150m","SPOT_890m",
                                   "SPOT_890m","SPOT_890m","SPOT_890m","SPOT_surface","SPOT_surface",
                                   "SPOT_surface","SPOT_surface","SPOT_surface","SPOT_surface",
                                   "SPOT_surface","SPOT_surface",
                                   "SPOT_surface","SPOT_surface",
                                   "SPOT_surface","SPOT_surface"),
                        "Replicate"=c("Rep1andRep2","Rep1andRep2","Rep1","Rep2","Rep1","Rep2",
                                      "Rep1andRep2","Rep1andRep2","Rep1","Rep2","Rep1","Rep2","19","20",
                                      "21","22","23","24","1","2","3","4","5","6","Rep1and2",
                                      "Rep3and4","Rep5and6","Rep1and2","Rep3and4","Rep5and6","Rep7",
                                      "7","8","9","10","11","12","13","14","!5","16","17","18"))


sample_merged <- data.frame("Files"=files) %>% 
  tidyr::separate(Files,sep="salmon_sample/",into=c("Stem","Name")) %>%
  tidyr::separate(Name,sep="_quant",into=c("SName","Extra")) %>% 
  dplyr::select(SName) %>%
  dplyr::left_join(sample_table,by=c("SName"="SRR")) %>%
  tidyr::unite("Sample_rep", Sample, Replicate, sep = "_", remove = FALSE)

names(files) <- sample_merged$Sample_rep

# header <- read.table("/../final-files/02-annotation_table/TaxonomicAndFunctionalAnnotations.csv", header = TRUE, nrow = 1,sep=",")

# tx2gene_2 <- data.frame(fread("/../final-files/02-annotation_table/TaxonomicAndFunctionalAnnotations.csv",skip=1,header=FALSE))
# 
# 
# tx2gene_2 <- data.frame(fread("/vortexfs1/omics/alexander/data/Hu-2022-ALOHA-SPOT/Hu_et_al_2022_eukrhythmic/final-files/02-annotation_table/TaxonomicAndFunctionalAnnotations.csv",skip=1,header=FALSE))

tx2gene_2 <- data.frame(fread("/../TaxonomicAndFunctionalAnnotations.csv", header = TRUE))

# ?setnames()
setnames(tx2gene_2, colnames(header))

ptm <- proc.time()

tx2gene_in=tx2gene_2 %>% dplyr::select(ShortSeqID,SequenceID)

txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene_in)

rm(tx2gene_2)
rm(tx2gene_in)

print(proc.time()-ptm)

library(tidyverse)

cat("\nPrep txi and merged sample table for DEseq input\n")

# Also make a sample table
metadata <- read_delim("complete-sample-list.txt", delim = ",")
sample_merged <- left_join(sample_merged, metadata, by = c("SName" = "RUN"))
write_delim(sample_merged, file = "sample_merged_txi.txt")

sample_merged <- read_delim("sample_merged_txi.txt")
rownames(sample_merged) <- sample_merged$Sample_rep
rownames(sample_merged) <- colnames(txi$counts)

save(txi, sample_merged, file = "/vortexfs1/scratch/sarahhu/tximport_metaT-ALOHA-CA.Rdata")

# Saves txi object, as this step requires a lot of memory to run.

## For subsetting
npsg_only
tmp_sample_merged <- sample_merged %>% 
  filter(Sample_rep %in% as.character(npsg_only$sample))
rownames(tmp_sample_merged) <- tmp_sample_merged$Sample_rep
rownames(tmp_sample_merged) <- colnames(txi_npsg$counts)