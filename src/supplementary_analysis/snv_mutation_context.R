library(tidyverse)
library(MutationalPatterns)
library(cowplot)

### Read in data 
df <- readr::read_delim("data/processed/mutational_signatures/input/output/SBS/mPPGL.SBS96.all", delim = "\t")
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

### Create germline category
sdh <- c("SDHB", "SDHA", "SDHC")
meta <- dplyr::mutate(meta, germline_cat = case_when(germline %in% sdh ~ "SDHx", TRUE ~ "Non-SDHx"))

### Create sample lists
samples_all <- meta$sample
samples_primary <- meta[meta$category == "Primary",]$sample
samples_met <- meta[meta$category == "Metastatic",]$sample
samples_pheo <- meta[meta$tumor_type == "PCC",]$sample
samples_para <- meta[meta$tumor_type == "PGL",]$sample
samples_sdh <- meta[meta$germline_cat == "SDHx",]$sample
samples_nosdh <- meta[meta$germline_cat == "Non-SDHx",]$sample

### Define plotting function
create_plot <- function(data, samples) {
  
  df_types <- data %>%
    dplyr::mutate(Mut = str_sub(MutationType, 3, 5)) %>%
    dplyr::select(-MutationType) %>%
    dplyr::group_by(Mut) %>%
    dplyr::summarise(across(everything(), sum)) %>%
    tibble::column_to_rownames(var = "Mut") %>% 
    dplyr::select(all_of(samples)) %>%
    t() %>% as.data.frame()
  
  plt <- plot_spectrum(df_types)
  plt <- plt + theme(axis.text.y = element_text(size=8), 
                     axis.title.y = element_text(size = 8),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6))
  return(plt)
  
}

### Create plots
plt1 <- create_plot(df, samples_all)
plt2 <- create_plot(df, samples_primary)
plt3 <- create_plot(df, samples_met)
plt4 <- create_plot(df, samples_pheo)
plt5 <- create_plot(df, samples_para)
plt6 <- create_plot(df, samples_sdh)
plt7 <- create_plot(df, samples_nosdh)

figure <- plot_grid(plt1 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt2 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt3 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt4 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt5 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt6 + theme(legend.position = "none", axis.title.y = element_blank()), 
                    plt7 + theme(legend.position = "none", axis.title.y = element_blank()), ncol=1)

pdf("results/figures/mutational_signatures/mutation_context.pdf", width = 7.5, height = 8)
print(figure)
dev.off()

# Plot 96 Mutation Context ------------------------------------------------

create_96_channel_plot <- function(data, samples) {
  df <- data %>% column_to_rownames(var = "MutationType") %>% 
    dplyr::select(all_of(samples)) %>%
    rowSums() %>% as.matrix()
  
  plt <- plot_96_profile(df, condensed = TRUE)
  plt <- plt + theme(strip.background.y = element_blank(), strip.text.y = element_blank())
  return(plt)
}

### Create plots
plt1 <- create_96_channel_plot(df, samples_all)
plt2 <- create_96_channel_plot(df, samples_primary)
plt3 <- create_96_channel_plot(df, samples_met)
plt4 <- create_96_channel_plot(df, samples_pheo)
plt5 <- create_96_channel_plot(df, samples_para)
plt6 <- create_96_channel_plot(df, samples_sdh)
plt7 <- create_96_channel_plot(df, samples_nosdh)

figure <- plot_grid(plt1 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt3 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt4 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt5 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt6 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), 
                    plt7 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank()), ncol=1)

pdf("results/figures/mutational_signatures/snv_mutation_96_channel_context.pdf", width = 6, height = 8.5)
print(figure)
dev.off()
