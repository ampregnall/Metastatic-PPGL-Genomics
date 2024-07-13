library(tidyverse)
library(ggpubr)
library(cowplot)

### Load data
df <- readxl::read_xlsx("metadata/Sequenza-Ploidy-Estimates.xlsx")
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
df <- df %>% dplyr::left_join(meta, by = "sample")

### Create SDH category
sdh <- c("SHDB", "SDHA", "SDHC")
df <- dplyr::mutate(df, germline_cat = case_when(germline %in% sdh ~ "SDHx", TRUE ~ "Non-SDHx"))

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
          axis.text.x = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.line.y = element_line())
}

# Create plots
plot1 <- create_plot(df, "category", "cellularity", "category", "Purity", "")
plot2 <- create_plot(df, "tumor_type", "cellularity", "tumor_type", "", "")
plot3 <- create_plot(df, "germline_cat", "cellularity", "germline_cat", "", "")

# Create overall plot
figure_a <- cowplot::plot_grid(plot1, plot2, plot3, nrow = 1)


# Ploidy Plots ------------------------------------------------------------

# Create the plots using the function
plot4 <- create_plot(df, "category", "ploidy", "category", "Ploidy", "")
plot5 <- create_plot(df, "tumor_type", "ploidy", "tumor_type", "", "")
plot6 <- create_plot(df, "germline_cat", "ploidy", "germline_cat", "", "")

# Create overall plot
figure_b <- cowplot::plot_grid(plot4, plot5, plot6, nrow = 1)


# Create Final Figure -----------------------------------------------------

figure <- cowplot::plot_grid(figure_a, figure_b, ncol = 1, align = "h", labels = "AUTO")

### Save results
pdf("results/figures/cnv_analysis/purity_ploidy_comparison.pdf", width = 16, height = 8.5)
print(figure)
dev.off()
