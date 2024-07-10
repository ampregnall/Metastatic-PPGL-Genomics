library(tidyverse)

### Load data
df <- readr::read_delim("data/processed/mutation_timing/ccf_estimates/mPPGL_cancer_cell_fractions.txt")

### Load IBD scores files for convenience --> contains tumor pairings
jsi <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "paired")

### Define empty vectors to store mutation information
metastasis_private_clonal <- vector(mode = "numeric", length = nrow(jsi))
primary_private_clonal <- vector(mode = "numeric", length = nrow(jsi))
metastasis_private_subclonal <- vector(mode = "numeric", length = nrow(jsi))
primary_private_subclonal <- vector(mode = "numeric", length = nrow(jsi))
shared_subclonal <- vector(mode = "numeric", length = nrow(jsi))
shared_clonal <- vector(mode = "numeric", length = nrow(jsi))

### Loop through sample pairs to calculate number of private/shared clonal/subclonal mutations
for (row in 1:nrow(jsi)) {
  ### Select mutations for each sample
  tmp1 <- df %>% dplyr::filter(Tumor.ID == jsi[row, "sample1"]$sample1)
  tmp2 <- df %>% dplyr::filter(Tumor.ID == jsi[row, "sample2"]$sample2)
  
  ### Count number of private clonal mutations in metastasis
  metastasis_private_clonal_priv <- tmp2 %>% dplyr::filter(CCF.95 == 1 & !(MutID %in% tmp1$MutID))
  metastasis_private_clonal[row] <- nrow(metastasis_private_clonal_priv)
  
  ### Count number of private clonal mutations in primary
  primary_private_clonal_priv <- tmp1 %>% dplyr::filter(CCF.95 == 1 & !(MutID %in% tmp2$MutID))
  primary_private_clonal[row] <- nrow(primary_private_clonal_priv)
  
  ### Count number of private subclonal mutations in primary
  metastasis_private_subclonal_priv <- tmp2 %>% dplyr::filter(CCF.95 < 1 & !(MutID %in% tmp1$MutID))
  metastasis_private_subclonal[row] <- nrow(metastasis_private_subclonal_priv)
  
  ### Count number of private clonal mutations in primary
  primary_private_subclonal_priv <- tmp1 %>% dplyr::filter(CCF.95 < 1 & !(MutID %in% tmp2$MutID))
  primary_private_subclonal[row] <- nrow(primary_private_subclonal_priv)
  
  ### Count number of shared subclonal mutation in paired samples
  shared_mutations <- intersect(tmp1$MutID, tmp2$MutID)
  shared_subclonal_priv <- tmp1 %>% dplyr::filter(CCF.95 < 1 & MutID %in% shared_mutations)
  shared_clonal_priv <- tmp1 %>% dplyr::filter(CCF.95 == 1 & MutID %in% shared_mutations)
  
  shared_subclonal[row] <- nrow(shared_subclonal_priv)
  shared_clonal[row] <- nrow(shared_clonal_priv)
}

### Add data to dataframe
jsi$primary_private_clonal <- primary_private_clonal
jsi$metastasis_private_clonal <- metastasis_private_clonal
jsi$metastasis_private_subclonal <- metastasis_private_subclonal
jsi$primary_private_subclonal <- primary_private_subclonal
jsi$shared_subclonal <- shared_subclonal
jsi$shared_clonal <- shared_clonal

### Calculate Jaccard similarity index
jsi <- jsi %>% dplyr::rowwise() %>%
  dplyr::mutate(JSI = shared_subclonal / (metastasis_private_clonal + primary_private_clonal + shared_subclonal))

### Save results
readr::write_delim(jsi, "data/processed/mutation_timing/jaccard_similarity/mPPGL_jsi_estimates.txt", delim = "\t")
