library(tidyverse)
library(ComplexHeatmap)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
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
mut_matrix = mut_matrix[, -1]
mut_matrix = t(as.matrix(mut_matrix))

### Define colors for Oncoplot
#col_cohort <- c("Derivation" = "#696DA1", "Validation"="#D745A1")
#col_category <- c("Primary" = "#2D3A7F", "Metastatic" = "lightgrey")
col_mutation = c("STOPGAIN" = "#A62C30", "MISSENSE" = "#9CC2C7", "FRAMESHIFT" = "#3075AC", "SPLICE" = "#A4C38E",
                 "SYNONYMOUS" = "#4E9858", "INTRON" = "#D98F8E", "INFRAME" = "#D4383C", "EXTRAGENIC" = "#E08244")
#col_tumor <- c("PCC" = "red", "PGL"="green")

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

oncoPrint(lof_matrix, alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, column_order = tumor_metadata$sample,
          top_annotation = HeatmapAnnotation(Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                             Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort, ),
                                             Type = tumor_metadata$category,
                                             Germline = tumor_metadata$germline))

### Create oncoPrint
plot <- oncoPrint(mut_matrix, alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, column_order = tumor_metadata$sample,
          top_annotation = HeatmapAnnotation(mutCount = anno_barplot(tumor_metadata$n_mutations, border = FALSE),
                                             Tumor = tumor_metadata$tumor_type,
                                             Cohort = tumor_metadata$cohort, 
                                             Type = tumor_metadata$category,
                                             Germline = tumor_metadata$germline))

### Save results
pdf("results/figures/driver_analysis/snv_oncoplot.pdf", width = 16, height = 8.5)
print(plot)
dev.off()
