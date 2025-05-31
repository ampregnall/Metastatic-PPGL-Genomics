## ---------------------------
## Script Name: misc_stats.R
## Description: Misc. analyses for manuscript
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-02-01
## Copyright (c) Andrew Pregnall, 2025
## ---------------------------

library(readr)
library(tidyverse)

# Test rates of DNA damage response and chromatin remodeling --------------

# Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

# Load driver mutations
lof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
drivers <- rbind(lof, gof)

drivers <- drivers %>%
  mutate(Variant.Type = 1)

# Get matrix driver of mutations in mPPGL cohort
mut_matrix <- drivers %>%
  distinct(Tumor.ID, Gene, .keep_all = TRUE) %>%
  select(Tumor.ID, Gene, Variant.Type) %>%
  pivot_wider(names_from = Gene, id_cols = Tumor.ID, values_from = Variant.Type, values_fill = 0)

# Some samples have no driver mutations; add back into matrix
mut_matrix <- left_join(samples, mut_matrix, by = c("sample" = "Tumor.ID"))
mut_matrix[is.na(mut_matrix)] <- 0


mut_matrix <- mut_matrix %>% 
  select(sample, KMT2D, KMT2C, ATRX, CREBBP, KMT2A, ATM, BRCA1, BRCA2, ATR) %>%
  mutate(mutation = case_when(rowSums(.[2:10]) > 0 ~ 1, TRUE ~ 0))

mut_matrix <- left_join(tumor_metadata, mut_matrix, by = c("sample" = "sample"))
mut_matrix <- mut_matrix %>% 
  mutate(germline_cat = case_when(germline %in% c("SDHB") ~ "SDHB", TRUE ~ "Non-SDHB"))

table(mut_matrix$germline_cat, mut_matrix$mutation)
fisher.test(table(mut_matrix$germline_cat, mut_matrix$mutation))

# Get rates of mutations in genes reported by Fishbein et al (2017) -------

tcga.genes <- c("HRAS", "NF1", "EPAS1", "RET", "CSDE1", "SETD2", "VHL", "FGFR1",
                "TP53", "BRAF", "ATRX", "ARNT", "IDH1")

drivers.fil <- drivers %>% filter(Gene %in% tcga.genes)

drivers.fil %>% count(Gene)
drivers.fil %>% distinct(Tumor.ID, Gene, .keep_all = TRUE) %>% count(Gene)

library(tidyverse)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

# Extract patient ID from the Tumor.ID and get unique values
tumor_metadata <- tumor_metadata %>%
  dplyr::rowwise() %>%
  dplyr::mutate(patient = str_split(sample, "-")[[1]][1]) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(patient, germline)

### Load driver information
lof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
snvs <- rbind(lof_snvs, gof_snvs)

### How many patients/tumors with mutations?
muts_all <- snvs %>% distinct(Normal.ID, Gene) %>%
  dplyr::filter(Gene %in% c("BRCA1", "BRCA2", "ATM", "ATR", "ATRX", "KMT2A", "KMT2C", "KMT2D")) %>%
  dplyr::distinct(Normal.ID) %>% dplyr::pull()

muts_chr <- snvs %>% distinct(Normal.ID, Gene) %>%
  dplyr::filter(Gene %in% c("ATRX", "KMT2A", "KMT2C", "KMT2D")) %>%
  dplyr::distinct(Normal.ID) %>% dplyr::pull()

muts_ddr <- snvs %>% distinct(Normal.ID, Gene) %>%
  dplyr::filter(Gene %in% c("BRCA1", "BRCA2", "ATM", "ATR")) %>%
  dplyr::distinct(Normal.ID) %>% dplyr::pull()

tumor_metadata <- tumor_metadata %>% dplyr::rowwise() %>%
  dplyr::mutate(muts_all = if_else(any(str_detect(muts_all, patient)), 1, 0),
                muts_chr = if_else(any(str_detect(muts_chr, patient)), 1, 0),
                muts_ddr = if_else(any(str_detect(muts_ddr, patient)), 1, 0),
                germ_cat = if_else(str_detect(germline, "SDH"), "SDH", "Non-SDH")) 

# Overall difference
table(tumor_metadata$muts_all, tumor_metadata$germ_cat)
fisher.test(table(tumor_metadata$muts_all, tumor_metadata$germ_cat))

# Difference in chromatin remodeling
table(tumor_metadata$muts_chr, tumor_metadata$germ_cat)
fisher.test(table(tumor_metadata$muts_chr, tumor_metadata$germ_cat))

# Difference in dna damage response
table(tumor_metadata$muts_ddr, tumor_metadata$germ_cat)
fisher.test(table(tumor_metadata$muts_ddr, tumor_metadata$germ_cat))

# Count of patients and tumors w/ mutations -------------------------------

### Number of tumors with mutations in given genes
snvs_tumor <- snvs %>% dplyr::distinct(Tumor.ID, Gene) %>%
  dplyr::group_by(Gene) %>%
  dplyr::count() %>%
  dplyr::mutate(percent = n / 48)

### Number of patients with mutations in given genes
snvs_patient <- snvs %>% dplyr::distinct(Normal.ID, Gene) %>%
  dplyr::group_by(Gene) %>%
  dplyr::count() %>%
  dplyr::mutate(percent = n / 27)
                                                                     