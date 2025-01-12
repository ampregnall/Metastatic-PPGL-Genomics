## ---------------------------
## Script Name: 04_snv_oncoplot.R
## Description: Generate Oncoplot of mPPGL driver mutations
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
##
## Date Created: 2025-01-11
## ---------------------------

# Load packages and define variables --------------------------------------

library(tidyverse)
library(ComplexHeatmap)

# Define variables
read_depth <- 8
allele_freq <- 0.01
panel_size <- 35.7

# Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

# Calculate tumor mutational burden ---------------------------------------

snvs <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.txt")

snvs <- snvs %>%
  filter(Tumor.AltDepth >= read_depth) %>%
  filter(gnomAD.MAX_AF == "." | as.numeric(gnomAD.MAX_AF) < allele_freq) %>%
  filter(EXON != ".") %>% # Select exonic mutations
  filter(!str_detect(Variant.Consequence, "synonymous_variant")) # Remove synonymous mutations

tmb <- snvs %>%
  group_by(Tumor.ID) %>%
  summarise(n_mutations = n()) %>%
  mutate(tmb = n_mutations / panel_size)

tumor_metadata <- left_join(tumor_metadata, tmb[, -2], by = c("sample" = "Tumor.ID"))

# Calculate number of driver mutations per sample  -------------------------

# Load driver information
lof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
drivers <- rbind(lof, gof)

### Calculate mutation burden of LOF mutations
driver_count <- drivers %>%
  group_by(Tumor.ID) %>%
  summarise(n_mutations = n())

# Merge information and calculate number of driver mutations
tumor_metadata <- left_join(tumor_metadata, driver_count, by = c("sample" = "Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] <- 0

# Prepare mutational signature data ---------------------------------------

# Plot signatures ---------------------------------------------------------

# Load SBS signatures
sigs <- read_delim("data/processed/mutational_signatures/SigProfilerExtractor_SBS/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt", delim = "\t")
sigs <- as.data.frame(sigs)
rownames(sigs) <- sigs[, 1]
sigs <- sigs[, -1]
sigs <- sigs[match(samples$sample, rownames(sigs)), ]
sigs <- t(apply(sigs, 1, function(x) x / sum(x))) # Normalize to percent

# Load indel signatures
ids <- readr::read_delim("data/processed/mutational_signatures/SigProfilerExtractor_ID/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities.txt", delim = "\t")
ids <- left_join(samples, ids, by = c("sample" = "Samples"))
ids <- as.data.frame(ids)
rownames(ids) <- ids[, 1]
ids <- ids[, -1]
ids <- ids[match(samples$sample, rownames(ids)), ]
ids <- t(apply(ids, 1, function(x) x / sum(x))) # Normalize to percent
ids[is.na(ids)] <- 0


# Create Matrix of LOF/GOF Mutations ------------------------------------------

# Create variant classification for Oncoplot
drivers <- drivers %>%
  mutate(Variant.Type = case_when(
    str_detect(Variant.Consequence, "stop_gained") ~ "STOPGAIN",
    str_detect(Variant.Consequence, "missense_variant") ~ "MISSENSE",
    str_detect(Variant.Consequence, "frameshift_variant") ~ "FRAMESHIFT",
    str_detect(Variant.Consequence, "splice") ~ "SPLICE"
  ))

# Get matrix driver of mutations in mPPGL cohort
mut_matrix <- drivers %>%
  distinct(Tumor.ID, Gene, .keep_all = TRUE) %>%
  select(Tumor.ID, Gene, Variant.Type) %>%
  pivot_wider(names_from = Gene, id_cols = Tumor.ID, values_from = Variant.Type)

# Some samples have no driver mutations; add back into matrix
mut_matrix <- left_join(samples, mut_matrix, by = c("sample" = "Tumor.ID"))

# Limit matrix to genes mutated in at least 6 samples
mut_matrix <- mut_matrix[colSums(is.na(mut_matrix)) <= 42]

# Transform into Oncoplot input
mut_matrix <- as.data.frame(mut_matrix)
mut_matrix[is.na(mut_matrix)] <- ""
rownames(mut_matrix) <- mut_matrix[, 1]
mut_matrix <- mut_matrix[match(samples$sample, rownames(mut_matrix)), ]
mut_matrix <- mut_matrix[, -1]
mut_matrix <- t(as.matrix(mut_matrix))

### Define colors for Oncoplot
category <- c(
  "Primary" = "#2D3A7F",
  "Metastatic" = "lightgrey"
)

mutation <- c(
  "STOPGAIN" = "#A62C30",
  "MISSENSE" = "#9CC2C7",
  "FRAMESHIFT" = "#3075AC",
  "SPLICE" = "#A4C38E"
)

tumor <- c(
  "PCC" = "#879eb3",
  "PGL" = "#ad823a"
)

germline <- c(
  "SDHA" = "#D3672D",
  "SDHB" = "#1B4367",
  "SDHC" = "#5F8528",
  "RET" = "#A9361E",
  "WT" = "#3F8B87"
)

alter_fun <- list(
  background = alter_graphic("rect", fill = "lightgrey"),
  STOPGAIN = alter_graphic("rect", fill = mutation["STOPGAIN"]),
  MISSENSE = alter_graphic("rect", fill = mutation["MISSENSE"]),
  FRAMESHIFT = alter_graphic("rect", fill = mutation["FRAMESHIFT"]),
  SPLICE = alter_graphic("rect", fill = mutation["SPLICE"])
)

# Set graphical parameters
ht_opt$COLUMN_ANNO_PADDING <- unit(4, "mm")

### Create oncoPrint
plot1 <- oncoPrint(mut_matrix,
  alter_fun = alter_fun, col = mutation, show_column_names = FALSE,
  column_order = tumor_metadata$sample, show_heatmap_legend = FALSE,
  top_annotation = HeatmapAnnotation(
    TMB = anno_barplot(tumor_metadata$tmb, height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))),
    Driver = anno_barplot(tumor_metadata$n_mutations, height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))),
    Germline = anno_simple(tumor_metadata$germline, col = germline),
    Tumor = anno_simple(x = tumor_metadata$tumor_type, col = tumor),
    Type = anno_simple(x = tumor_metadata$category, col = category),
    show_legend = FALSE
  ),
  bottom_annotation = HeatmapAnnotation(
    SBS = anno_barplot(sigs,
      gp = gpar(fill = c(
        "#948958",
        "#798259",
        "#febf70",
        "#a17721",
        "#9f6a4d",
        "#a95035",
        "#96211e"
      )),
      height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))
    ),
    Indels = anno_barplot(ids,
      gp = gpar(fill = c(
        "#5a8f9b",
        "#d0d7da",
        "#4c6a74",
        "#f2e6d7"
      )),
      height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))
    ),
    gap = unit(4, "mm")
  )
)


### Save results
pdf("results/figures/driver_analysis/snv_oncoplot.pdf", width = 15, height = 12)
print(plot1)
dev.off()


# Oncoplot Sensitivity Analysis -------------------------------------------

# Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

# Calculate score for number of callers that support a mutation
drivers <- drivers %>%
  mutate(LANCET = case_when(LANCET == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(MUTECT2 = case_when(MUTECT2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(STRELKA2 = case_when(STRELKA2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(VARDICT = case_when(VARDICT == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(VARSCAN2 = case_when(VARSCAN2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  rowwise() %>%
  mutate(SCORE = LANCET + MUTECT2 + STRELKA2 + VARDICT + VARSCAN2) %>%
  filter(!(Variant.Class != "SNV" & SCORE < 2))

# Calculate mutation burden of LOF mutations
driver_count <- drivers %>%
  group_by(Tumor.ID) %>%
  summarise(n_mutations = n())

# Merge information and calculate number of driver mutations
tumor_metadata <- left_join(tumor_metadata, driver_count, by = c("sample" = "Tumor.ID"))
tumor_metadata <- left_join(tumor_metadata, tmb[, -2], by = c("sample" = "Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] <- 0

# Get matrix of variants in cohort
mut_matrix <- drivers %>%
  distinct(Tumor.ID, Gene, .keep_all = TRUE) %>%
  select(Tumor.ID, Gene, Variant.Type) %>%
  pivot_wider(names_from = Gene, id_cols = Tumor.ID, values_from = Variant.Type)

# Some samples have no mutations; add back into matrix
mut_matrix <- left_join(samples, mut_matrix, by = c("sample" = "Tumor.ID"))

### Limit matrix to genes mutated in at least 6 samples
mut_matrix <- mut_matrix[colSums(is.na(mut_matrix)) <= 42]

### Transform into Oncoplot input
mut_matrix <- as.data.frame(mut_matrix)
mut_matrix[is.na(mut_matrix)] <- ""
rownames(mut_matrix) <- mut_matrix[, 1]
mut_matrix <- mut_matrix[match(samples$sample, rownames(mut_matrix)), ]
mut_matrix <- mut_matrix[, -1]
mut_matrix <- t(as.matrix(mut_matrix))

### Create oncoPrint
plot2 <- oncoPrint(mut_matrix,
  alter_fun = alter_fun, col = mutation, show_column_names = FALSE,
  column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
  top_annotation = HeatmapAnnotation(
    TMB = anno_barplot(tumor_metadata$tmb, height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))),
    Driver = anno_barplot(tumor_metadata$n_mutations, height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))),
    Germline = anno_simple(tumor_metadata$germline, col = germline),
    Tumor = anno_simple(x = tumor_metadata$tumor_type, col = tumor),
    Type = anno_simple(x = tumor_metadata$category, col = category),
    show_legend = TRUE
  )
)

### Save results
pdf("results/figures/driver_analysis/snv_oncoplot_sensitivity_analysis.pdf", width = 15, height = 12)
print(plot2)
dev.off()
