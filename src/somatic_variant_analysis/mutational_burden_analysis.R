library(tidyverse)
library(ggpubr)
library(cowplot)

# Read in somatic variant data
df1 <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.txt")
df2 <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.txt")
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

### Calculate number of LOF mutations per sample
df1_summary <- df1 %>% dplyr::group_by(Tumor.ID) %>%
  dplyr::summarise(n_lof = n())

### Filter for GOF mutations and calculate mutations per sample
df2_summary <- df1 %>% dplyr::group_by(Tumor.ID) %>%
  dplyr::summarise(n_gof = n())

### Join data frames into summary for plotting
muts_per_sample <- dplyr::left_join(tumor_metadata, df1_summary, by=c("sample" = "Tumor.ID"))
muts_per_sample <- muts_per_sample %>% dplyr::left_join(df2_summary, by=c("sample" = "Tumor.ID")) 
muts_per_sample$n_total <- muts_per_sample$n_gof + muts_per_sample$n_lof
muts_per_sample[is.na(muts_per_sample)] = 0

### Create SDH category
sdh <- c("SDHB", "SDHA", "SDHC")
muts_per_sample <- dplyr::mutate(muts_per_sample, germline_cat = case_when(germline %in% sdh ~ "SDHx", TRUE ~ "Non-SDHx"))

# Purity Plots ------------------------------------------------------------

# Define a function to create box plots
create_plot <- function(data, x_var, y_var, fill_var, y_label, x_label) {
  ggboxplot(data = data, x = x_var, y = y_var, fill = fill_var, 
            palette = c('#083D77', "#F4D35E"), ylab = y_label, 
            xlab = x_label, add = "jitter", 
            add.params = list(color = "#071E22", size = 2)) +
    stat_compare_means(method = "wilcox.test", label.x = 1.5, size = 7, 
                       aes(label = paste0("p = ", after_stat(p.format)))) +
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

# Create plots
plot1 <- create_plot(muts_per_sample, "category", "n_total", "category", "Count", "")
plot2 <- create_plot(muts_per_sample, "tumor_type", "n_total", "tumor_type", "", "")
plot3 <- create_plot(muts_per_sample, "germline_cat", "n_total", "germline_cat", "", "")

plot2 <- plot2 + theme(axis.text.y = element_blank())
plot3 <- plot3 + theme(axis.text.y = element_blank())


wilcox.test(n_total ~ category, data = muts_per_sample)
wilcox.test(n_total ~ tumor_type, data = muts_per_sample)
wilcox.test(n_total ~ germline_cat, data = muts_per_sample)

# Create overall plot
figure <- cowplot::plot_grid(plot1, plot2, plot3, nrow = 1)

### Save results
pdf("results/figures/driver_analysis/mutational_burden_comparison.pdf", width = 7.5, height = 5)
print(figure)
dev.off()