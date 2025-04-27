## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-10-06
## Copyright (c) Andrew Pregnall, 2024
## ---------------------------

library(tidyverse)

### Load data
df <- readr::read_delim("data/processed/cnvs/annotSV/mPPGL.annotSV.data.combined.txt")
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
genes <- readr::read_delim("metadata/gene_with_protein_product.txt")

### Get number of metastatic and primary tumors
n_metastatic <- sum(meta$category == "Metastatic")
n_primary <- sum(meta$category == "Primary")

### extract metastatic tumor samples; get number of samples
mets <- meta %>% 
  dplyr::filter(category == "Metastatic") %>%
  dplyr::pull(sample)

primaries <- meta %>% 
  dplyr::filter(category == "Primary") %>%
  dplyr::pull(sample)

### Merge sample information onto cnv data
df <- df %>% dplyr::select(TumorID, Gene, Segment.Type) %>%
  dplyr::mutate(Segment.Type.Score = case_when(
    Segment.Type == "amp" ~ 2,
    Segment.Type == "gain" ~ 1,
    Segment.Type == "neutral" ~ 0,
    Segment.Type == "loss" ~ -1,
    Segment.Type == "del" ~ -2
  ))

### Calculate enrichment data for primaries
primary_muts <- df %>% 
  dplyr::filter(TumorID %in% primaries & Gene %in% genes$symbol) %>%
  dplyr::distinct(TumorID, Gene, .keep_all = TRUE) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(score = sum(Segment.Type.Score) / n_primary) %>%
  dplyr::arrange(desc(score))

### Save results
readr::write_delim(
  x = as.data.frame(primary_muts), delim = "\t", col_names = FALSE,
  file = "data/processed/enrichment_analysis/primary_cnv_mutations.rnk"
)
    
### Calculate enrichment data for metastases
met_muts <- df %>% 
  dplyr::filter(TumorID %in% mets & Gene %in% genes$symbol) %>%
  dplyr::distinct(TumorID, Gene, .keep_all = TRUE) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(score = sum(Segment.Type.Score) / n_metastatic) %>%
  dplyr::arrange(desc(score))

### Save results
readr::write_delim(
  x = as.data.frame(met_muts), delim = "\t", col_names = FALSE,
  file = "data/processed/enrichment_analysis/metastatic_cnv_mutations.rnk"
)
