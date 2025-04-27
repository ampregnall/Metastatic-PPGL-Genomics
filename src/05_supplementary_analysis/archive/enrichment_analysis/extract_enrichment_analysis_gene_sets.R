## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-09-27
## ---------------------------

### Load required libraries
library(tidyverse)

### Load data and merge
df <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.txt", delim = "\t")
df2 <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.txt", delim = "\t")
all <- rbind(df, df2)

### Load metadata
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
genes <- readr::read_delim("metadata/gene_with_protein_product.txt")

### extract primary tumor samples; get number of samples
primaries <- meta %>% dplyr::filter(category == "Primary") %>%
  dplyr::pull(sample)

n_primaries <- length(primaries)

### extract metastatic tumor samples; get number of samples
mets <- meta %>% dplyr::filter(category == "Metastatic") %>%
  dplyr::pull(sample)

n_mets <- length(mets)

met_muts <- all %>% dplyr::filter(Tumor.ID %in% mets & Gene %in% genes$symbol) %>%
  dplyr::distinct(Tumor.ID, Gene) %>%
  dplyr::count(Gene) %>%
  dplyr::mutate(n = n / n_mets) %>%
  dplyr::arrange(desc(n))

primary_muts <- all %>% dplyr::filter(Tumor.ID %in% primaries & Gene %in% genes$symbol) %>%
  dplyr::distinct(Tumor.ID, Gene) %>%
  dplyr::count(Gene) %>%
  dplyr::mutate(n = n / n_primaries) %>%
  dplyr::arrange(desc(n))


readr::write_delim(x = as.data.frame(met_muts), delim = "\t", col_names = FALSE,
                   file = "data/processed/enrichment_analysis/metastasis_mutations.rnk")

readr::write_delim(x = as.data.frame(primary_muts), delim = "\t", col_names = FALSE,
                   file = "data/processed/enrichment_analysis/primary_mutations.rnk")

