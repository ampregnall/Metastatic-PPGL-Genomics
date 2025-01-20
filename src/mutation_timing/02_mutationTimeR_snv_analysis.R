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
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Gene", y = "Number of Mutations", fill = "Mutation Type") +
  ggtitle("Proportion of Timing Estimates for Driver Mutations") +
  theme(legend.position = "bottom")

pdf("results/figures/mutational_timing/mutation_timing_driver_mutations.pdf", width = 10, height = 5)
print(plt4)
dev.off()
