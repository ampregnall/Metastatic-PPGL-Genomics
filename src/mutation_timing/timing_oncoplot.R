library(tidyverse)
library(ComplexHeatmap)
library(circlize)

### Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
tumor_metadata <- tumor_metadata %>% dplyr::arrange(order)
samples <- tumor_metadata %>% dplyr::select(sample)

# Clonal/subclonal information --------------------------------------------
df <- readr::read_delim("data/processed/mutation_timing/ccf_estimates/mPPGL_cancer_cell_fractions.txt")

### Classify mutations as clonal/subclonal based on CCF
df <- df %>% dplyr::mutate(CLS = case_when(CCF.95 == 1 ~ "Clonal", TRUE ~ "Subclonal"))

### Count number of clonal/subclonal mutations by sample
CLS_sample_summary <- df %>% dplyr::group_by(Tumor.ID, CLS) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = CLS, values_from = count) %>%
  tidyr::replace_na(list(Subclonal = 0)) %>%
  dplyr::mutate(Prop.Clonal = Clonal / (Subclonal + Clonal),
                Prop.Subclonal = Subclonal / (Subclonal + Clonal)) %>%
  dplyr::select(Tumor.ID, Prop.Clonal, Prop.Subclonal)

CLS_sample_summary <- as.data.frame(CLS_sample_summary)
rownames(CLS_sample_summary) = CLS_sample_summary[, 1]
CLS_sample_summary <- CLS_sample_summary[match(samples$sample, rownames(CLS_sample_summary)), ]
CLS_sample_summary = CLS_sample_summary[, -1]

### Load copy number information
df$Chr <- as.character(df$Chr)
cytos <- readr::read_delim("metadata/GrCh38.arm.coordinates.txt")
cytos <- cytos %>% dplyr::rename(Cyto.Start = Start, Cyto.End = End)

arms <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
          "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q",
          "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q",
          "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q")

### MERGE CYTOBAND INFORMATION AND SELECT APPROPRIATE ARM
df <- df %>% dplyr::full_join(cytos, by = c("Chr" = "Chrom"), relationship = "many-to-many") %>% 
  dplyr::filter(SV.Start >= Cyto.Start & SV.Start < Cyto.End) %>%
  dplyr::mutate(Arm = factor(Arm, levels = arms))


df_timing <- df %>% dplyr::group_by(Tumor.ID, Arm, timing) %>%
  dplyr::summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = timing, values_from = count, values_fill = 0) %>%
  dplyr::rowwise() %>% dplyr::mutate(timing_prop = late / (early + late)) %>%
  dplyr::select(Tumor.ID, Arm, timing_prop) %>%
  tidyr::pivot_wider(names_from = Arm, values_from = timing_prop)

df_timing <- as.data.frame(df_timing)
rownames(df_timing) = df_timing[, 1]
df_timing <- df_timing[match(samples$sample, rownames(df_timing)), ]
df_timing = df_timing[, -1]

### Add missing chromosomes
df_timing$`13p` <- NA
df_timing$`14p` <- NA
df_timing$`15p` <- NA
df_timing$`22p` <- NA

col_fun = colorRamp2(c(0, 0.5, 1), c("#3c7b4a", "white", "#6b3e8a"))
### Define colors for Oncoplot
col_cohort <- c("Derivation" = "#696DA1", 
                "Validation"=  "#D745A1")

col_category <- c("Primary" = "#2D3A7F", 
                  "Metastatic" = "lightgrey")

col_tumor <- c("PCC" = "#879eb3", 
               "PGL"="#ad823a")

col_germline <- c("SDHA" = "#D3672D", 
                  "SDHB" = "#1B4367", 
                  "SDHC" = "#5F8528", 
                  "RET" =  "#A9361E", 
                  "WT" =   "#3F8B87")


### Create plot
plot <- Heatmap(as.matrix(df_timing),
  col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE, column_order = arms,
  gap = 1, show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = TRUE,
  left_annotation = rowAnnotation(Tumor = anno_simple(x = tumor_metadata$tumor_type, col = col_tumor),
                                  Cohort = anno_simple(x = tumor_metadata$cohort, col = col_cohort),
                                  Type = anno_simple(x = tumor_metadata$category, col = col_category),
                                  Germline = anno_simple(tumor_metadata$germline, col = col_germline), show_legend = FALSE),
  right_annotation = rowAnnotation(Clonality = anno_barplot(CLS_sample_summary,
    gp = gpar(fill = c("lightblue", "#D36492"), border = FALSE), width = unit(1, "in")
  ))
)

### Save results
pdf("results/figures/mutational_timing/timing_oncoplot.pdf", width = 15, height = 10)
print(plot)
dev.off()



