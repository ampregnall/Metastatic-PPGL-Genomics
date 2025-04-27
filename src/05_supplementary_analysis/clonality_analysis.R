library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

# Driver Mutation Clonality Analysis ------------------------------------------------------

### Load clonality data
df <- readr::read_delim("data/processed/mutation_timing/ccf_estimates/mPPGL_cancer_cell_fractions.txt", delim = "\t")

### Load driver mutations
lof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
snvs <- rbind(lof_snvs, gof_snvs)

### Load metadata and extract whether samples are pheo or para
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
pheos <- dplyr::filter(meta, tumor_type == "PCC")$sample
paras <- dplyr::filter(meta, tumor_type == "PGL")$sample

### Load list of paired samples and limit to primary-metastasis pairs
tumor_pairs <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "paired")

### Add mutation ID to SNV data then left join clonality data for driver mutations
snvs <- snvs %>% dplyr::mutate(MutID = stringr::str_c(Chr, Start, Gene, REF, ALT, sep = "-")) %>%
  dplyr::left_join(df, by = c("Tumor.ID", "MutID"))

### Get number of driver mutations lacking CCF estimate; would annoTSV files reduce?
sum(is.na(snvs$CCF)) # 208 w/ Sequenza; 90 w/ annotSV! (82 muts on X chromosome)

### Calculate proportion of shared, primary-private, and metastasis-private clonal drivers
snvs <- snvs %>% dplyr::filter(!is.na(CCF.95)) %>%
  dplyr::mutate(CLS = case_when(CCF.95 == 1 ~ "Clonal", TRUE ~ "Subclonal"))

### Define empty vectors to store mutation information
metastasis_private_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
primary_private_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
metastasis_private_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
primary_private_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
shared_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
shared_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))

### Loop through sample pairs to calculate number of private/shared clonal/subclonal mutations
for (row in 1:nrow(tumor_pairs)) {
  ### Select mutations for each sample
  tmp1 <- snvs %>% dplyr::filter(Tumor.ID == tumor_pairs[row, "sample1"]$sample1)
  tmp2 <- snvs %>% dplyr::filter(Tumor.ID == tumor_pairs[row, "sample2"]$sample2)
  
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
tumor_pairs$primary_private_clonal <- primary_private_clonal
tumor_pairs$metastasis_private_clonal <- metastasis_private_clonal
tumor_pairs$metastasis_private_subclonal <- metastasis_private_subclonal
tumor_pairs$primary_private_subclonal <- primary_private_subclonal
tumor_pairs$shared_subclonal <- shared_subclonal
tumor_pairs$shared_clonal <- shared_clonal

### Select primary-met pairs only
primary_met_pairs <- tumor_pairs %>% dplyr::filter(type == "primary_metastasis")
primary_met_pairs_pheo <- primary_met_pairs %>% dplyr::filter(sample1 %in% pheos)
primary_met_pairs_para <- primary_met_pairs %>% dplyr::filter(sample1 %in% paras)

tbl <- primary_met_pairs %>% dplyr::select(primary_private_clonal, primary_private_subclonal, 
                                     metastasis_private_clonal,metastasis_private_subclonal) %>%
  pivot_longer(cols = everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(n = sum(value))

tbl_pheo <- primary_met_pairs_pheo %>% dplyr::select(primary_private_clonal, primary_private_subclonal, 
                                           metastasis_private_clonal,metastasis_private_subclonal) %>%
  pivot_longer(cols = everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(n = sum(value))

tbl_para <- primary_met_pairs_para %>% dplyr::select(primary_private_clonal, primary_private_subclonal, 
                                                     metastasis_private_clonal,metastasis_private_subclonal) %>%
  pivot_longer(cols = everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(n = sum(value))

### Define helper function to calculate proportion of primary private clonal mutations, 
### metastasis private clonal mutations, and shared clonal mutations
calculate_clonal_proportions <- function(df, label) {
  tmp <- df %>% dplyr::select(primary_private_clonal, metastasis_private_clonal, shared_clonal) %>%
    colSums() %>% as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "type") %>%
    dplyr::mutate(proportion = value / sum(value))
  tmp$group <- label
  return(tmp)
} 

### Calculate proportions as above
clonal_prop_all <- calculate_clonal_proportions(primary_met_pairs, "All")
clonal_prop_pheo <- calculate_clonal_proportions(primary_met_pairs_pheo, "Pheo.")
clonal_prop_para <- calculate_clonal_proportions(primary_met_pairs_para, "Para.") 

### Bind data and label
clonal_prop <- rbind(clonal_prop_all, clonal_prop_para, clonal_prop_pheo)
clonal_prop$facet_label <- "Clonal drivers"

### Define helper function to calculate proportion of primary private subclonal mutations, 
### metastasis private subclonal mutations, and shared subclonal mutations
calculate_subclonal_proportions <- function(df, label) {
  tmp <- df %>% dplyr::select(primary_private_subclonal, metastasis_private_subclonal, shared_subclonal) %>%
    colSums() %>% as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "type") %>%
    dplyr::mutate(proportion = value / sum(value))
  tmp$group <- label
  return(tmp)
}

### Calculate proportions as above
subclonal_prop_all <- calculate_subclonal_proportions(primary_met_pairs, "All")
subclonal_prop_pheo <- calculate_subclonal_proportions(primary_met_pairs_pheo, "Pheo.")
subclonal_prop_para <- calculate_subclonal_proportions(primary_met_pairs_para, "Para.") 

### Bind data and label
subclonal_prop <- rbind(subclonal_prop_all, subclonal_prop_para, subclonal_prop_pheo)
subclonal_prop$facet_label <- "Subclonal drivers"

### Create plotting data
props <- rbind(clonal_prop, subclonal_prop)

### Create categorical group for mutation type
metastasis_private <- c("metastasis_private_clonal", "metastasis_private_subclonal")
primary_private <- c("primary_private_clonal", "primary_private_subclonal")
shared <- c("shared_clonal", "shared_subclonal")

props <- props %>% dplyr::mutate(type_grouped = case_when(type %in% metastasis_private ~ "Metastasis-private",
                                                          type %in% primary_private ~ "Primary-private",
                                                          type %in% shared ~ "Shared"))

### Create plot
plot <- ggplot(props, aes(x = group, y = proportion, fill = type_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ facet_label, nrow=1) +
  labs(title = "",
       x = "", y = "", fill = "Type") +
  scale_fill_manual(values = c("#DF6050", "#4BA789", "#62A1CA"),
                    labels = c("Metastasis-private", "Primary-private", "Shared")) +
  theme_minimal() + theme(plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 16),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 16), 
                          legend.position = "none",
                          strip.text = element_blank())

### Save results
pdf("results/figures/clonality_analysis/driver_clonality.pdf", width = 3.75, height = 1.5)
print(plot)
dev.off()

# All Mutation Clonality Analysis -----------------------------------------

df <- readr::read_delim("data/processed/mutation_timing/jaccard_similarity/mPPGL_jsi_estimates.txt", delim = "\t")

### Select primary-met pairs only
primary_met_pairs <- df %>% dplyr::filter(type == "primary_metastasis")
primary_met_pairs_pheo <- df %>% dplyr::filter(sample1 %in% pheos)
primary_met_pairs_para <- df %>% dplyr::filter(sample1 %in% paras)

### Calculate proportions as above
clonal_prop_all <- calculate_clonal_proportions(primary_met_pairs, "All")
clonal_prop_pheo <- calculate_clonal_proportions(primary_met_pairs_pheo, "Pheo.")
clonal_prop_para <- calculate_clonal_proportions(primary_met_pairs_para, "Para.") 

### Bind data and label
clonal_prop <- rbind(clonal_prop_all, clonal_prop_para, clonal_prop_pheo)
clonal_prop$facet_label <- "Clonal mutations"

### Calculate proportions as above
subclonal_prop_all <- calculate_subclonal_proportions(primary_met_pairs, "All")
subclonal_prop_pheo <- calculate_subclonal_proportions(primary_met_pairs_pheo, "Pheo.")
subclonal_prop_para <- calculate_subclonal_proportions(primary_met_pairs_para, "Para.") 

### Bind data and label
subclonal_prop <- rbind(subclonal_prop_all, subclonal_prop_para, subclonal_prop_pheo)
subclonal_prop$facet_label <- "Subclonal mutations"

### Create plotting data
props <- rbind(clonal_prop, subclonal_prop)

props <- props %>% dplyr::mutate(type_grouped = case_when(type %in% metastasis_private ~ "Metastasis-private",
                                                          type %in% primary_private ~ "Primary-private",
                                                          type %in% shared ~ "Shared"))

### Create plot
plot <- ggplot(props, aes(x = group, y = proportion, fill = type_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ facet_label, nrow=1) +
  labs(title = "",
       x = "", y = "", fill = "Type") +
  scale_fill_manual(values = c("#DF6050", "#4BA789", "#62A1CA"),
                    labels = c("Metastasis-private", "Primary-private", "Shared")) +
  theme_minimal() + theme(plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 16),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 16), 
                          legend.position = "none",
                          strip.text = element_blank())

### Save results
pdf("results/figures/clonality_analysis/mutation_clonality.pdf", width = 3.75, height = 1.5)
print(plot)
dev.off()

# Jaccard Statistics ------------------------------------------------------

### Load data
df <- readr::read_delim("data/processed/mutation_timing/jaccard_similarity/mPPGL_jsi_estimates.txt", delim = "\t")

### Add annotation for monoclonal versus polyclonal seeding
df <- df %>% dplyr::filter(type == "primary_metastasis") %>% 
  dplyr::mutate(Seeding = case_when(JSI >= 0.3 ~ "Polyclonal", TRUE ~ "Monoclonal"))

### Calculate proportion of metastases that exhibited monoclonal vs polyclonal seeding
seeding_summary <- df %>% dplyr::group_by(Seeding) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(proportion = count / sum(count))

### Calculate interquartile range of JSI by monoclonal vs polyclonal seeding status
aggregate(JSI~Seeding, data=df, summary)

### Compare number of certain mutations in monoclonal vs polyclonal seeding states
wilcox.test(primary_private_clonal ~ Seeding, data=df)
aggregate(primary_private_clonal ~ Seeding, data=df, summary)

wilcox.test(metastasis_private_clonal ~ Seeding, data=df)
aggregate(metastasis_private_clonal ~ Seeding, data=df, summary)

wilcox.test(shared_subclonal ~ Seeding, data=df)
aggregate(shared_subclonal ~ Seeding, data=df, summary)

df_p <- dplyr::filter(df, Seeding == "Polyclonal")
df_m <- dplyr::filter(df, Seeding == "Monoclonal")

### Create plot
p1 <- ggplot(df_m, aes(x = type, y=JSI)) + 
  geom_beeswarm(size = 3, cex = 3, color = "#2E86AB") +
  geom_beeswarm(data = df_p, aes(type, y = JSI), size = 3, cex = 3, color = "#564138") +
  xlab("") + ylim(0, 0.6) + theme_minimal() + 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_blank(), 
        axis.text.x = element_blank())

### Save results
pdf("results/figures/clonality_analysis/jaccard_beeswarm.pdf", width = 6, height = 4.5)
print(p1)
dev.off()

