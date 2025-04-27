# Oncoplot Sensitivity Analysis -------------------------------------------

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
  filter(!str_detect(Variant.Consequence, "synonymous_variant")) %>%
  mutate(LANCET = case_when(LANCET == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(MUTECT2 = case_when(MUTECT2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(STRELKA2 = case_when(STRELKA2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(VARDICT = case_when(VARDICT == "PASS" ~ 1, TRUE ~ 0)) %>%
  mutate(VARSCAN2 = case_when(VARSCAN2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  rowwise() %>%
  mutate(SCORE = LANCET + MUTECT2 + STRELKA2 + VARDICT + VARSCAN2) %>%
  filter(!(Variant.Class != "SNV" & SCORE < 2))# Remove synonymous mutations

tmb <- snvs %>%
  group_by(Tumor.ID) %>%
  summarise(n_mutations = n()) %>%
  mutate(tmb = n_mutations / panel_size)

tumor_metadata <- left_join(tumor_metadata, tmb[, -2], by = c("sample" = "Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] <- 0


# Process driver mutations ------------------------------------------------

# Load driver information
lof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof <- read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
drivers <- rbind(lof, gof)

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

drivers <- drivers %>%
  mutate(Variant.Type = case_when(
    str_detect(Variant.Consequence, "stop_gained") ~ "STOPGAIN",
    str_detect(Variant.Consequence, "missense_variant") ~ "MISSENSE",
    str_detect(Variant.Consequence, "frameshift_variant") ~ "FRAMESHIFT",
    str_detect(Variant.Consequence, "splice") ~ "SPLICE"
  ))

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
  "PGL" = "#542E71"
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
ht_opt$ROW_ANNO_PADDING <- unit(4, "mm")

### Create oncoPrint
plot2 <- oncoPrint(mut_matrix,
                   alter_fun = alter_fun, col = mutation, show_column_names = FALSE,
                   column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
                   top_annotation = HeatmapAnnotation(
                     TMB = anno_barplot(tumor_metadata$tmb, height = unit(1, "in"), border = FALSE, axis_param = list(gp = gpar(fontsize = 12))),
                     Germline = anno_simple(tumor_metadata$germline, col = germline),
                     Tumor = anno_simple(x = tumor_metadata$tumor_type, col = tumor),
                     Type = anno_simple(x = tumor_metadata$category, col = category),
                     show_legend = TRUE
                   ))

### Save results
pdf("results/figures/driver_analysis/snv_oncoplot_sensitivity_analysis.pdf", width = 15, height = 12)
print(plot2)
dev.off()