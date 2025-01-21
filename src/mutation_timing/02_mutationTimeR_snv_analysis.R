## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
##
## Date Created: 2025-01-13
## ---------------------------

# Load packages
library(tidyverse)

# Define variables
col <- c("#4DAF4A", "#984EA3", "#377EB8", "#E41A1C", "#999999")
read_depth <- 8
allele_freq <- 0.01
clin_var_filter <- c("Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic")
drivers <- c(
  "KMT2A", "KMT2D", "TRIP11", "AKAP9", "ZFHX3", "ATR", "TPR", "ATM", "GRIN2A", "POLR2A", "BRCA1",
  "MUC16", "ASXL1", "SETD2", "TET2", "KMT2C", "ATRX", "CHD2", "LRP1B", "POLQ", "UBR5", "CSMD3", "MTOR",
  "STAG1", "CREBBP", "ROS1", "SPEN", "FAT3", "BRCA2", "MYH9", "NSD1", "ASPM", "MUC4", "KAT6B", "FAT1",
  "CIC"
)

# Define variables
clonal <- c("clonal [late]", "clonal [NA]", "clonal [early]")
subclonal <- c("subclonal")

# Load data -------
files <- list.files(
  path = "data/processed/mutation_timing/MutationTimeR/variant_timing",
  pattern = "*.csv", full.names = FALSE, recursive = FALSE
)

df <- map(files, function(file) {
  read_csv(paste0("data/processed/mutation_timing/MutationTimeR/variant_timing/", file),
    col_types = cols(
      ...1 = col_double(),
      Tumor.ID = col_character(),
      Normal.ID = col_character(),
      Chr = col_character(),
      Start = col_double(),
      REF = col_character(),
      ALT = col_character(),
      FILTER = col_character(),
      Gene = col_character(),
      Gene.Accession = col_character(),
      Variant.Class = col_character(),
      Variant.Type = col_character(),
      Variant.Consequence = col_character(),
      HGVSc = col_character(),
      HGVSp = col_character(),
      Feature.Type = col_character(),
      Feature.Accession = col_character(),
      Bio.type = col_character(),
      Existing.variation = col_character(),
      EXON = col_character(),
      INTRON = col_character(),
      STRAND = col_double(),
      cDNA.position = col_character(),
      CDS.position = col_character(),
      Protein.position = col_character(),
      Amino.acids = col_character(),
      Codons = col_character(),
      SpliceAI.DS_AG = col_character(),
      SpliceAI.DS_AL = col_character(),
      SpliceAI.DS_DG = col_character(),
      SpliceAI.DS_DL = col_character(),
      REVEL = col_character(),
      gnomAD.AF = col_character(),
      gnomAD.AFR = col_character(),
      gnomAD.AMR = col_character(),
      gnomAD.ASJ = col_character(),
      gnomAD.EAS = col_character(),
      gnomAD.FIN = col_character(),
      gnomAD.NFE = col_character(),
      gnomAD.OTH = col_character(),
      gnomAD.SAS = col_character(),
      gnomAD.MAX_AF = col_character(),
      gnomAD.MAX_POPS = col_character(),
      LOFTEE.lof = col_character(),
      LOFTEE.filter = col_character(),
      LOFTEE.flags = col_character(),
      LOFTEE.info = col_character(),
      ClinVar = col_character(),
      ClinVar.SIG = col_character(),
      ClinVar.REVSTAT = col_character(),
      ClinVar.DN = col_character(),
      Tumor.Zyg = col_character(),
      Tumor.Depth = col_double(),
      Tumor.AltDepth = col_double(),
      Tumor.AltFrac = col_double(),
      Normal.Zyg = col_character(),
      Normal.Depth = col_double(),
      Normal.AltDepth = col_double(),
      Normal.AltFrac = col_double(),
      LANCET = col_character(),
      MUTECT2 = col_character(),
      STRELKA2 = col_character(),
      VARDICT = col_character(),
      VARSCAN2 = col_character(),
      MutCN = col_double(),
      MutDeltaCN = col_double(),
      MajCN = col_double(),
      MinCN = col_double(),
      MajDerCN = col_double(),
      MinDerCN = col_double(),
      CNF = col_double(),
      CNID = col_double(),
      pMutCN = col_double(),
      pGain = col_double(),
      pSingle = col_double(),
      pSub = col_double(),
      pMutCNTail = col_double(),
      CLS = col_character(),
      CN = col_double(),
      chromosome = col_character(),
      start.pos = col_double(),
      end.pos = col_double(),
      Bf = col_double(),
      N.BAF = col_double(),
      sd.BAF = col_double(),
      depth.ratio = col_double(),
      N.ratio = col_double(),
      sd.ratio = col_double(),
      CNt = col_double(),
      A = col_double(),
      B = col_double(),
      LPP = col_double()
    )
  )
}) %>% bind_rows()

# Overall Distribution for SBS and Indels ---------------------------------

# Get proportion of different mutations in each sample
df.prop <- df %>%
  group_by(Tumor.ID) %>%
  count(CLS) %>%
  group_by(Tumor.ID) %>%
  mutate(prop = n / sum(n))

# Plot distribution of timing estimates for all mutations
plt1 <- ggplot(df.prop, aes(x = Tumor.ID, y = prop, fill = CLS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Sample ID", y = "Proportion of Mutations", fill = "Mutation Type") +
  ggtitle("Proportion of Timing Estimates in Each Sample") +
  theme(legend.position = "bottom")

pdf("results/figures/mutational_timing/mutation_timing_proportions_all.pdf", width = 10, height = 5)
print(plt1)
dev.off()

# Subset data by mutation type
df.snv <- df %>% filter(Variant.Class == "SNV")
df.id <- df %>% filter(Variant.Class != "SNV")

# Get proportion of different mutations in each sample for SNVs
df.snv.prop <- df.snv %>%
  group_by(Tumor.ID) %>%
  count(CLS) %>%
  group_by(Tumor.ID) %>%
  mutate(prop = n / sum(n))

# Plot distribution of timing estimates for SNVs
plt2 <- ggplot(df.snv.prop, aes(x = Tumor.ID, y = prop, fill = CLS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Sample ID", y = "Proportion of Mutations", fill = "Mutation Type") +
  ggtitle("Proportion of Timing Estimates in Each Sample (SNVs)") +
  theme(legend.position = "bottom")

pdf("results/figures/mutational_timing/mutation_timing_proportions_snv.pdf", width = 10, height = 5)
print(plt2)
dev.off()

# Get proportion of different mutations in each sample for indels
df.id.prop <- df.id %>%
  group_by(Tumor.ID) %>%
  count(CLS) %>%
  group_by(Tumor.ID) %>%
  mutate(prop = n / sum(n))

# Plot distribution of timing estimates for indels
plt3 <- ggplot(df.id.prop, aes(x = Tumor.ID, y = prop, fill = CLS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Sample ID", y = "Proportion of Mutations", fill = "Mutation Type") +
  ggtitle("Proportion of Timing Estimates in Each Sample (Indels)") +
  theme(legend.position = "bottom")

pdf("results/figures/mutational_timing/mutation_timing_proportions_indel.pdf", width = 10, height = 5)
print(plt3)
dev.off()

# Driver Mutations --------------------------------------------------------

# Load CGC mutations
cgc <- read_delim("metadata/cancer_census_genes_all_v98.tsv", delim = "\t")

# Get driver mutations
df.drivers <- df %>%
  filter(Tumor.AltDepth >= read_depth) %>%
  filter(!str_detect(tolower(ClinVar.SIG), "benign")) %>%
  filter(gnomAD.MAX_AF == "." | as.numeric(gnomAD.MAX_AF) < allele_freq) %>%
  mutate(
    mutation.consequence =
      case_when(
        str_detect(Variant.Consequence, "frameshift_variant|stop_gained") ~ "LOF",
        Variant.Consequence %in% c("splice_acceptor_variant", "splice_donor_variant") ~ "LOF",
        str_detect(Variant.Consequence, "missense_variant") &
          REVEL != "." & as.numeric(REVEL) > 0.5 ~ "LOF",
        str_detect(Variant.Consequence, "missense_variant") &
          (ClinVar.SIG %in% clin_var_filter | grepl("COSV", Existing.variation)) ~ "GOF",
        TRUE ~ "Unknown"
      )
  ) %>%
  filter(mutation.consequence != "Unknown") %>%
  filter(Gene %in% cgc$`Gene Symbol`)

#
df.drivers.prop <- df.drivers %>%
  filter(Gene %in% drivers) %>%
  group_by(Gene) %>%
  count(CLS) %>%
  group_by(Gene) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  inner_join(df.drivers %>% group_by(Gene, CLS) %>% count(), by = "Gene")

# Plot distribution of timing estimates for driver mutations
plt4 <- ggplot(df.drivers.prop, aes(x = reorder(Gene, -total), y = n, fill = CLS)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, size = 8),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "Count", fill = "") +
  ggtitle("Proportion of Timing Estimates for Driver Mutations") +
  theme(legend.position = "none")

pdf("results/figures/mutational_timing/mutation_timing_driver_mutations.pdf", width = 7.5, height = 1.5)
print(plt4)
dev.off()

### Load metadata and extract whether samples are pheo or para
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
pheos <- dplyr::filter(meta, tumor_type == "PCC")$sample
paras <- dplyr::filter(meta, tumor_type == "PGL")$sample

df.drivers <- df.drivers %>% mutate(MutID = str_c(Chr, Start, Gene, REF, ALT, sep = "-"))

### Load list of paired samples and limit to primary-metastasis pairs
tumor_pairs <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "paired")

### Define empty vectors to store mutation information
metastasis_private_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
primary_private_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
metastasis_private_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
primary_private_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
shared_subclonal <- vector(mode = "numeric", length = nrow(tumor_pairs))
shared_clonal <- vector(mode = "numeric", length = nrow(tumor_pairs))

### Loop through sample pairs to calculate number of private/shared clonal/subclonal mutations
for (row in 1:nrow(tumor_pairs)) {
  # Select mutations for each sample
  tmp1 <- df.drivers %>% dplyr::filter(Tumor.ID == tumor_pairs[row, "sample1"]$sample1)
  tmp2 <- df.drivers %>% dplyr::filter(Tumor.ID == tumor_pairs[row, "sample2"]$sample2)
  
  ### Count number of private clonal mutations in metastasis
  metastasis_private_clonal_priv <- tmp2 %>% dplyr::filter(CLS %in% clonal & !(MutID %in% tmp1$MutID))
  metastasis_private_clonal[row] <- nrow(metastasis_private_clonal_priv)
  
  ### Count number of private clonal mutations in primary
  primary_private_clonal_priv <- tmp1 %>% dplyr::filter(CLS %in% clonal & !(MutID %in% tmp2$MutID))
  primary_private_clonal[row] <- nrow(primary_private_clonal_priv)
  
  ### Count number of private subclonal mutations in primary
  metastasis_private_subclonal_priv <- tmp2 %>% dplyr::filter(CLS %in% subclonal & !(MutID %in% tmp1$MutID))
  metastasis_private_subclonal[row] <- nrow(metastasis_private_subclonal_priv)
  
  ### Count number of private subclonal mutations in primary
  primary_private_subclonal_priv <- tmp1 %>% dplyr::filter(CLS %in% subclonal & !(MutID %in% tmp2$MutID))
  primary_private_subclonal[row] <- nrow(primary_private_subclonal_priv)
  
  ### Count number of shared subclonal mutation in paired samples
  shared_mutations <- intersect(tmp1$MutID, tmp2$MutID)
  shared_subclonal_priv <- tmp1 %>% dplyr::filter(CLS %in% subclonal & MutID %in% shared_mutations)
  shared_clonal_priv <- tmp1 %>% dplyr::filter(CLS %in% clonal & MutID %in% shared_mutations)
  
  shared_subclonal[row] <- nrow(shared_subclonal_priv)
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
plt5 <- ggplot(props, aes(x = group, y = proportion, fill = type_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ facet_label, nrow=1) +
  labs(title = "",
       x = "", y = "", fill = "Type") +
  scale_fill_manual(values = c("#DF6050", "#4BA789", "#62A1CA"),
                    labels = c("Metastasis-private", "Primary-private", "Shared")) +
  theme_minimal() + theme(plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(size = 16),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 16), 
                          legend.position = "none",
                          strip.text = element_text(size = 8))

### Save results
pdf("results/figures/mutational_timing/driver_clonality.pdf", width = 5.4, height = 1.7)
print(plt5)
dev.off()