## ---------------------------
## Script Name: cnv_oncoplot.R
## Description: Create CNV oncoplot for mPPGL cohort
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-09-27
## ---------------------------

library(tidyverse)
library(ComplexHeatmap)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

### Load copy number information
cnvs_arm_level <- readr::read_delim("data/processed/cnvs/arm_level_changes/mPPGL_arm_level_summary.txt")
cnvs <- readr::read_delim("data/processed/cnvs/annotSV/mPPGL.annotSV.data.combined.txt")


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
  dplyr::mutate(Change = case_when(Arm.Gain == 1 ~ "Gain", Arm.Loss == 1 ~ "Loss", LOH == 1 ~ "Neutral LOH", TRUE ~ "")) %>%
  dplyr::select(Sample, arm, Change) %>%
  tidyr::pivot_wider(names_from = arm, values_from = Change)

### Transform into Oncoplot input
cnv_matrix <- as.data.frame(cnv_matrix)
rownames(cnv_matrix) = cnv_matrix[, 1]
cnv_matrix <- cnv_matrix[match(samples$sample, rownames(cnv_matrix)), ]
cnv_matrix = cnv_matrix[, -1]
cnv_matrix <- cnv_matrix[tumor_metadata$sample,,drop=FALSE]

### Define colors for Oncoplot
col_category <- c("Primary" = "#2D3A7F", 
                  "Metastatic" = "lightgrey")

col_tumor <- c("PCC" = "#879eb3", 
               "PGL"="#542E71")

col_germline <- c("SDHA" = "#D3672D", 
                  "SDHB" = "#1B4367", 
                  "SDHC" = "#5F8528", 
                  "RET" =  "#A9361E", 
                  "WT" =   "#3F8B87")

### Define colors for Oncoplot
col_mutation = c("Loss" = "#01016f", "Gain" = "#d8031c", "Neutral LOH" = "#b4d5e4")

alter_fun = list(
  background = alter_graphic("rect", fill = "lightgrey"),   
  Gain = alter_graphic("rect", fill = col_mutation["Gain"]),
  Loss = alter_graphic("rect", fill = col_mutation["Loss"]),
  'Neutral LOH' = alter_graphic("rect", fill = col_mutation["Neutral LOH"]))

ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")

### Create oncoPrint
plot <- oncoPrint(cnv_matrix, 
          alter_fun = alter_fun, col = col_mutation, show_column_names = TRUE, show_heatmap_legend = FALSE, 
          column_order = arms, row_order = tumor_metadata$sample, show_pct = FALSE, show_row_names = FALSE,
          left_annotation = rowAnnotation(Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                          Type = anno_simple(x = tumor_metadata$category, col = col_category, ),
                                          Germline = anno_simple(tumor_metadata$germline, col = col_germline), show_legend = FALSE),
          right_annotation = rowAnnotation(Count = anno_barplot(tumor_metadata$count, border = FALSE, width = unit(1, "in"), axis_param = list(gp = gpar(fontsize = 12)))))


### Save results
pdf("results/figures/cnv_analysis/cnv_oncoplot.pdf", width = 15, height = 10)
print(plot)
dev.off()