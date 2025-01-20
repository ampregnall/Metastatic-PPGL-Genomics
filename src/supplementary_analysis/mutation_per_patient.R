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
