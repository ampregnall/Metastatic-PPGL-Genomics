library(tidyverse)
library(ggpubr)
library(cowplot)

### Load data
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
cnvs <- readr::read_delim("data/processed/cnvs/annotSV/mPPGL.annotSV.data.combined.txt")

### Calculate number of copy number events per sample
cnv_count <- cnvs %>% dplyr::distinct(TumorID, AnnotSV.ID) %>% 
  dplyr::group_by(TumorID) %>% dplyr::summarise(count = n())

### Add purity, ploidy, and cnv count to data
meta <- dplyr::left_join(meta, cnv_count, by = c("sample"="TumorID"))
  
### Create SDH category
sdh <- c("SDHB")
meta <- dplyr::mutate(meta, germline_cat = case_when(germline %in% sdh ~ "SDHB", TRUE ~ "Non-SDHB"))

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

# Create plots
plot1 <- create_plot(meta, "category", "cellularity", "category", "Purity", "")
plot2 <- create_plot(meta, "tumor_type", "cellularity", "tumor_type", "", "")
plot3 <- create_plot(meta, "germline_cat", "cellularity", "germline_cat", "", "")

wilcox.test(cellularity ~ category, data = meta)
aggregate(meta$cellularity, list(meta$category), FUN = median)

wilcox.test(cellularity ~ tumor_type, data = meta)
aggregate(meta$cellularity, list(meta$tumor_type), FUN = median)


wilcox.test(cellularity ~ germline_cat, data = meta)


# Ploidy Plots ------------------------------------------------------------

# Create the plots using the function
plot4 <- create_plot(meta, "category", "ploidy", "category", "Ploidy", "")
plot5 <- create_plot(meta, "tumor_type", "ploidy", "tumor_type", "", "")
plot6 <- create_plot(meta, "germline_cat", "ploidy", "germline_cat", "", "")

### Perform stats
wilcox.test(ploidy ~ category, data = meta) # p = 0.1663
wilcox.test(ploidy ~ tumor_type, data = meta) # p = 0.7045
wilcox.test(ploidy ~ germline_cat, data = meta) # p = 0.3484

# Event comparison --------------------------------------------------------

### Perform stat testing before removing outlier for plotting
wilcox.test(count ~ category, data = meta) # p = 0.768
wilcox.test(count ~ tumor_type, data = meta) # p = 0.3412
wilcox.test(count ~ germline_cat, data = meta) # p = 0.7937

### Remove outlier for plot
meta <- meta %>% dplyr::filter(sample != "PP220-DZ1A")

# Create the plots using the function
plot7 <- create_plot(meta, "category", "count", "category", "Copy number events", "")
plot8 <- create_plot(meta, "tumor_type", "count", "tumor_type", "", "")
plot9 <- create_plot(meta, "germline_cat", "count", "germline_cat", "", "")

# Create Final Figure -----------------------------------------------------

figure <- cowplot::plot_grid(plot1, plot2, plot3,
                             plot4, plot5, plot6, 
                             plot7, plot8, plot9,
                             ncol = 3, nrow = 3, align = "hv")

### Save results
pdf("results/figures/cnv_analysis/purity_ploidy_event_comparison.pdf", width = 7.5, height = 7.5)
print(figure)
dev.off()
