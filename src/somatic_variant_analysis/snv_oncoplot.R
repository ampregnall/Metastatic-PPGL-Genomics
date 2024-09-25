library(tidyverse)
library(ComplexHeatmap)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

### Load driver information
lof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt")
gof_snvs <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt")
snvs <- rbind(lof_snvs, gof_snvs)

# Calculate Mutation Burden -----------------------------------------------

### Calculate mutation burden of LOF mutations
mutation_count <- snvs %>% dplyr::group_by(Tumor.ID) %>%
  dplyr::summarise(n_mutations = n())

### Merge information and calculate number of driver mutations
tumor_metadata <- dplyr::left_join(tumor_metadata, mutation_count, by=c("sample"="Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] <- 0

mutation_plot <- ggplot(tumor_metadata, aes(x=order, y=n_mutations)) + 
  geom_bar(stat = "identity") + scale_y_log10() +
  xlab("") + ylab("Count") + theme_minimal() + 
  theme(strip.background = element_blank(),
      strip.text = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),,
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      legend.text = element_text(size = 16),
      text = element_text(size = 16))

pdf("results/figures/driver_analysis/mutation_burden.pdf", width = 5, height = 0.75)
print(mutation_plot)
dev.off()

# Create Matrix of LOF/GOF Mutations ------------------------------------------

### Create variant classification for Oncoplot
snvs <- snvs %>% dplyr::mutate(Variant.Type = case_when(str_detect(Variant.Consequence, "stop_gained") ~ "STOPGAIN",
                                                                str_detect(Variant.Consequence, "missense_variant") ~ "MISSENSE",
                                                                str_detect(Variant.Consequence, "frameshift_variant") ~ "FRAMESHIFT",
                                                                str_detect(Variant.Consequence, "splice") ~ "SPLICE",
                                                                str_detect(Variant.Consequence, "synonymous_variant") ~ "SYNONYMOUS",
                                                                str_detect(Variant.Consequence, "inframe") ~ "INFRAME",
                                                                str_detect(Variant.Consequence, "UTR") ~ "UTR",
                                                                str_detect(Variant.Consequence, "stream") ~ "EXTRAGENIC",
                                                                str_detect(Variant.Consequence, "intron") ~ "INTRON",
                                                                TRUE ~ "OTHER"))

### Get matrix of GOF variants in cohort
mut_matrix <- snvs %>% dplyr::distinct(Tumor.ID, Gene, .keep_all = TRUE) %>% 
  dplyr::select(Tumor.ID, Gene, Variant.Type) %>%
  pivot_wider(names_from = Gene, id_cols = Tumor.ID, values_from = Variant.Type) 

### Some samples have no GOF mutations; add back into matrix
mut_matrix <- dplyr::left_join(samples, mut_matrix, by=c("sample"="Tumor.ID"))

### Limit matrix to genes mutated in at least 6 samples
mut_matrix <- mut_matrix[colSums(is.na(mut_matrix)) <= 42]

### Transform into Oncoplot input
mut_matrix <- as.data.frame(mut_matrix)
mut_matrix[is.na(mut_matrix)] = ""
rownames(mut_matrix) = mut_matrix[, 1]
mut_matrix <- mut_matrix[match(samples$sample, rownames(mut_matrix)), ]
mut_matrix = mut_matrix[, -1]
mut_matrix = t(as.matrix(mut_matrix))

### Define colors for Oncoplot
col_cohort <- c("Derivation" = "#696DA1", 
                "Validation"=  "#D745A1")

col_category <- c("Primary" = "#2D3A7F", "Metastatic" = "lightgrey")

col_mutation = c("STOPGAIN" = "#A62C30", "MISSENSE" = "#9CC2C7", "FRAMESHIFT" = "#3075AC", "SPLICE" = "#A4C38E",
                 "SYNONYMOUS" = "#4E9858", "INTRON" = "#D98F8E", "INFRAME" = "#D4383C", "EXTRAGENIC" = "#E08244")

col_tumor <- c("PCC" = "#879eb3", "PGL"="#ad823a")

col_germline <- c("SDHA" = "#D3672D", 
                  "SDHB" = "#1B4367", 
                  "SDHC" = "#5F8528", 
                  "RET" =  "#A9361E", 
                  "WT" =   "#3F8B87")

alter_fun = list(
  background = alter_graphic("rect", fill = "lightgrey"),   
  STOPGAIN = alter_graphic("rect", fill = col_mutation["STOPGAIN"]),
  MISSENSE = alter_graphic("rect", fill = col_mutation["MISSENSE"]),
  FRAMESHIFT = alter_graphic("rect", fill = col_mutation["FRAMESHIFT"]),
  SPLICE = alter_graphic("rect", fill = col_mutation["SPLICE"]),
  SYNONYMOUS = alter_graphic("rect", fill = col_mutation["SYNONYMOUS"]),
  INTRON = alter_graphic("rect", fill = col_mutation["INTRON"]),
  INFRAME = alter_graphic("rect", fill = col_mutation["INFRAME"]),
  EXTRAGENIC = alter_graphic("rect", fill = col_mutation["EXTRAGENIC"]))

### Create oncoPrint
plot1 <- oncoPrint(mut_matrix, alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, 
          column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
          top_annotation = HeatmapAnnotation(Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                             Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort),
                                             Type = anno_simple(x = tumor_metadata$category, col = col_category),
                                             Germline = anno_simple(tumor_metadata$germline, col = col_germline), 
                                             show_legend = FALSE))


### Create oncoPrint
plot2 <- oncoPrint(mut_matrix, alter_fun = alter_fun, col = col_mutation, 
                   show_column_names = FALSE, show_pct = FALSE, show_row_names = FALSE,
                   column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
                   top_annotation = HeatmapAnnotation(Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                                      Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort),
                                                      Type = anno_simple(x = tumor_metadata$category, col = col_category),
                                                      Germline = anno_simple(tumor_metadata$germline, col = col_germline), 
                                                      show_legend = FALSE))


### Save results
pdf("results/figures/driver_analysis/snv_oncoplot.pdf", width = 15, height = 10.5)
print(plot1)
print(plot2)
dev.off()

# Plot signatures ---------------------------------------------------------

### Load sigs
sigs <- readr::read_delim("results/mutational_signatures/SigProfilerExtractor/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt", delim = "\t")

sigs <- sigs %>% 
  reshape2::melt(id=c("Samples")) %>%
  dplyr::left_join(tumor_metadata, by=c("Samples"="sample"))

sbs_plot <- ggplot(sigs, aes(x=order, y=value, fill=variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#febf70", "#a17721", "#9f6a4d", "#a95035", "#96211e")) +
  xlab("") + ylab("Percent") + theme_minimal() + theme(strip.background = element_blank(),
                                                       strip.text = element_blank(),
                                                       panel.border = element_blank(),
                                                       panel.grid = element_blank(),
                                                       axis.text.x = element_blank(),
                                                       axis.text.y = element_blank(),,
                                                       axis.title.x = element_text(size = 16),
                                                       axis.title.y = element_blank(),
                                                       legend.text = element_text(size = 16),
                                                       legend.position = "none",
                                                       text = element_text(size = 16))

pdf("results/figures/driver_analysis/sbs_signatures.pdf", width = 5, height = 0.75)
print(sbs_plot)
dev.off()

sigs <- readr::read_delim("data/processed/mutational_signatures/SigProfilerExtractor_ID/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Activities/COSMIC_ID83_Activities.txt", delim = "\t")

sigs <- sigs %>% 
  reshape2::melt(id=c("Samples")) %>%
  dplyr::left_join(tumor_metadata, by=c("Samples"="sample"))

indel_plot <- ggplot(sigs, aes(x=order, y=value, fill=variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#5a8f9b", "#d0d7da", "#4c6a74", "#f2e6d7")) +
  xlab("") + ylab("Percent") + theme_minimal() + theme(strip.background = element_blank(),
                                                       strip.text = element_blank(),
                                                       panel.border = element_blank(),
                                                       panel.grid = element_blank(),
                                                       axis.text.x = element_blank(),
                                                       axis.text.y = element_blank(),,
                                                       axis.title.x = element_text(size = 16),
                                                       axis.title.y = element_blank(),
                                                       legend.text = element_text(size = 16),
                                                       legend.position = "none",
                                                       text = element_text(size = 16))

pdf("results/figures/driver_analysis/indel_signatures.pdf", width = 5, height = 0.75)
print(indel_plot)
dev.off()

# Oncoplot Sensitivity Analysis -------------------------------------------

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

snvs <- snvs %>% dplyr::mutate(LANCET = case_when(LANCET == "PASS" ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(MUTECT2 = case_when(MUTECT2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(STRELKA2 = case_when(STRELKA2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(VARDICT = case_when(VARDICT == "PASS" ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(VARSCAN2 = case_when(VARSCAN2 == "PASS" ~ 1, TRUE ~ 0)) %>%
  dplyr::rowwise() %>% dplyr::mutate(SCORE = LANCET + MUTECT2 + STRELKA2 + VARDICT + VARSCAN2) %>%
  dplyr::filter(!(Variant.Class != "SNV" & SCORE < 2))

### Calculate mutation burden of LOF mutations
mutation_count <- snvs %>% dplyr::group_by(Tumor.ID) %>%
  dplyr::summarise(n_mutations = n())

### Merge information and calculate number of driver mutations
tumor_metadata <- dplyr::left_join(tumor_metadata, mutation_count, by=c("sample"="Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] <- 0

### Get matrix of GOF variants in cohort
mut_matrix <- snvs %>% dplyr::distinct(Tumor.ID, Gene, .keep_all = TRUE) %>% 
  dplyr::select(Tumor.ID, Gene, Variant.Type) %>%
  pivot_wider(names_from = Gene, id_cols = Tumor.ID, values_from = Variant.Type) 

### Some samples have no GOF mutations; add back into matrix
mut_matrix <- dplyr::left_join(samples, mut_matrix, by=c("sample"="Tumor.ID"))

### Limit matrix to genes mutated in at least 6 samples
mut_matrix <- mut_matrix[colSums(is.na(mut_matrix)) <= 42]

### Transform into Oncoplot input
mut_matrix <- as.data.frame(mut_matrix)
mut_matrix[is.na(mut_matrix)] = ""
rownames(mut_matrix) = mut_matrix[, 1]
mut_matrix <- mut_matrix[match(samples$sample, rownames(mut_matrix)), ]
mut_matrix = mut_matrix[, -1]
mut_matrix = t(as.matrix(mut_matrix))

### Create oncoPrint
plot1 <- oncoPrint(mut_matrix, alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, 
          column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
          top_annotation = HeatmapAnnotation(Burden = anno_barplot(tumor_metadata$n_mutations, height = unit(1, "in"), border = FALSE),
                                             Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                             Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort),
                                             Type = anno_simple(x = tumor_metadata$category, col = col_category),
                                             Germline = anno_simple(tumor_metadata$germline, col = col_germline), 
                                             show_legend = FALSE))

### Create oncoPrint
plot2 <- oncoPrint(mut_matrix, alter_fun = alter_fun, col = col_mutation, 
                   show_column_names = FALSE, show_pct = FALSE, show_row_names = FALSE,
                   column_order = tumor_metadata$sample, show_heatmap_legend = FALSE, 
                   top_annotation = HeatmapAnnotation(Burden = anno_barplot(tumor_metadata$n_mutations, height = unit(1, "in"), border = FALSE),
                                                      Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                                      Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort),
                                                      Type = anno_simple(x = tumor_metadata$category, col = col_category),
                                                      Germline = anno_simple(tumor_metadata$germline, col = col_germline), 
                                                      show_legend = FALSE))

### Save results
pdf("results/figures/driver_analysis/snv_oncoplot_sensitivity_analysis.pdf", width = 15, height = 10.5)
print(plot1)
print(plot2)
dev.off()
