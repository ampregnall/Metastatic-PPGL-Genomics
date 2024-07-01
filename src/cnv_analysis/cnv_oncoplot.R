library(tidyverse)
library(ComplexHeatmap)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
samples <- tumor_metadata %>% dplyr::select(sample)

### Load copy number information
cnvs_arm_level <- readr::read_delim("data/processed/cnvs/PPGL.arm_level_changes.combined.txt")
cnvs <- readr::read_delim("data/processed/cnvs/PPGL.annotSV.data.combined.csv")


arms <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
          "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q",
          "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q",
          "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q", "Xp", "Xq")

# Calculate Mutation Burden -----------------------------------------------

cnv_count <- cnvs %>% dplyr::distinct(TumorID, AnnotSV.ID) %>% 
  dplyr::group_by(TumorID) %>% dplyr::summarise(count = n())

### Merge information and calculate number of driver mutations
tumor_metadata <- dplyr::left_join(tumor_metadata, cnv_count, by=c("sample"="TumorID"))
tumor_metadata[is.na(tumor_metadata)] <- 0

# Create Matrix of LOF/GOF Mutations ------------------------------------------

### Create variant classification for Oncoplot
cnv_matrix <- cnvs_arm_level %>% 
  dplyr::mutate(Change = case_when(Arm.Gain == 1 ~ "Gain", Arm.Loss == 1 ~ "Loss", TRUE ~ "")) %>%
  dplyr::select(Sample, Arm, Change) %>%
  tidyr::pivot_wider(names_from = Arm, values_from = Change)

### Transform into Oncoplot input
cnv_matrix <- as.data.frame(cnv_matrix)
rownames(cnv_matrix) = cnv_matrix[, 1]
cnv_matrix = cnv_matrix[, -1]
cnv_matrix_t = t(as.matrix(cnv_matrix))

test2 <- cnv_matrix[tumor_metadata$sample,,drop=FALSE]

### Define colors for Oncoplot
#col_cohort <- c("Derivation" = "#696DA1", "Validation"="#D745A1")
#col_category <- c("Primary" = "#2D3A7F", "Metastatic" = "lightgrey")
col_mutation = c("Loss" = "#01016f", "Gain" = "#d8031c")
#col_tumor <- c("PCC" = "red", "PGL"="green")

alter_fun = list(
  background = alter_graphic("rect", fill = "lightgrey"),   
  Gain = alter_graphic("rect", fill = col_mutation["Gain"]),
  Loss = alter_graphic("rect", fill = col_mutation["Loss"]))

### Create oncoPrint
plot <- oncoPrint(test2, alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, 
          column_order = arms, row_order = tumor_metadata$sample, show_pct = FALSE, show_row_names = FALSE,
          left_annotation = rowAnnotation(Sample = anno_text(tumor_metadata$sample),
                                          Tumor = tumor_metadata$tumor_type,
                                          Cohort = tumor_metadata$cohort,
                                          Type = tumor_metadata$category,
                                          Germline = tumor_metadata$germline),
          right_annotation = rowAnnotation(Count = anno_barplot(tumor_metadata$count, border = FALSE)))

### Save results
pdf("results/figures/cnv_analysis/cnv_oncoplot.pdf", width = 16, height = 8.5)
print(plot)
dev.off()



Bottom_ha = HeatmapAnnotation(Response = clin$Response, 
                              "TP53 Status" = clin$TP53_status,
                              col = list(Response = structure(c("purple4", "grey"), names = c("Responders", "NonResponders")),
                                         "TP53 Status" = structure(c("firebrick", "grey"), names = c("Mutant","Wild type"))),
                              border = TRUE)
