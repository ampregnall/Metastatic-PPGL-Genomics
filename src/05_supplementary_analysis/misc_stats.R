## ---------------------------
## Script Name: misc_stats.R
## Description: PURPOSE
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

                                                                     