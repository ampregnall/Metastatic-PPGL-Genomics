## ---------------------------
## Script Name: mutational_burden_analysis.R
## Description: Compares TMB and number of loss/gain of function variants in cohort
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-23
## ---------------------------

# Load packages
library(tidyverse)
library(ggpubr)
library(cowplot)

# Define variables
panel_size <- 35.7
read_depth <- 8
allele_freq <- 0.01
sdh <- c("SDHB")

# Load sample information
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

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

### Calculate number of LOF mutations per sample
drivers_summary <- drivers %>%
  group_by(Tumor.ID) %>%
  summarise(n_muts = n())

### Join data frames into summary for plotting
tumor_metadata <- left_join(tumor_metadata, drivers_summary, by = c("sample" = "Tumor.ID"))
tumor_metadata[is.na(tumor_metadata)] = 0

### Create SDH category
tumor_metadata <- mutate(tumor_metadata, germline_cat = case_when(germline %in% sdh ~ "SDHB", TRUE ~ "Non-SDHB"))

# Purity Plots ------------------------------------------------------------

# Define a function to create box plots
create_plot <- function(data, x_var, y_var, fill_var, y_label, x_label) {
  ggviolin(data = data, x = x_var, y = y_var, fill = fill_var, 
            palette = c('#083D77', "#F4D35E"), ylab = y_label, 
            xlab = x_label, add = "jitter", 
            add.params = list(color = "#071E22", size = 2)) +
    theme(legend.position = "none", 
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_text(size = 8, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.line.y = element_line())
}


# Plotting and statistical testing for tumor mutation burden --------------

# Create plots
plot1 <- create_plot(tumor_metadata, "category", "tmb", "category", "TMB", "")
plot2 <- create_plot(tumor_metadata, "tumor_type", "tmb", "tumor_type", "", "")
plot3 <- create_plot(tumor_metadata, "germline_cat", "tmb", "germline_cat", "", "")

# Edit legends
plot2 <- plot2 + theme(axis.text.y = element_blank())
plot3 <- plot3 + theme(axis.text.y = element_blank())

# Statistical testing
wilcox.test(tmb ~ category, data = tumor_metadata)
aggregate(x = tumor_metadata$tmb, by = list(tumor_metadata$category), FUN = median)

wilcox.test(tmb ~ tumor_type, data = tumor_metadata)
aggregate(x = tumor_metadata$tmb, by = list(tumor_metadata$tumor_type), FUN = median)

wilcox.test(tmb ~ germline_cat, data = tumor_metadata)
aggregate(x = tumor_metadata$tmb, by = list(tumor_metadata$germline_cat), FUN = median)


# Plotting and statistical testing for n drivers --------------------------

# Create plots
plot4 <- create_plot(tumor_metadata, "category", "n_muts", "category", "Drivers (Count)", "")
plot5 <- create_plot(tumor_metadata, "tumor_type", "n_muts", "tumor_type", "", "")
plot6 <- create_plot(tumor_metadata, "germline_cat", "n_muts", "germline_cat", "", "")

# Edit legends
plot5 <- plot5 + theme(axis.text.y = element_blank())
plot6 <- plot6 + theme(axis.text.y = element_blank())

# Statistical testing
wilcox.test(n_muts ~ category, data = tumor_metadata)
aggregate(x = tumor_metadata$n_muts, by = list(tumor_metadata$category), FUN = median)

wilcox.test(n_muts ~ tumor_type, data = tumor_metadata)
aggregate(x = tumor_metadata$n_muts, by = list(tumor_metadata$tumor_type), FUN = median)

wilcox.test(n_muts ~ germline_cat, data = tumor_metadata)
aggregate(x = tumor_metadata$n_muts, by = list(tumor_metadata$germline_cat), FUN = median)


# Create overall plot
figure <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 2)

### Save results
pdf("results/figures/driver_analysis/mutational_burden_comparison.pdf", width = 7.5, height = 7.5)
print(figure)
dev.off()
