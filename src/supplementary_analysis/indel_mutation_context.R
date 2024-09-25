library(tidyverse)
library(MutationalPatterns)
library(cowplot)

### Read in data 
df <- readr::read_delim("data/processed/mutational_signatures/input/output/ID/mPPGL.ID83.all", delim = "\t")
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

## Read in data from MutationPatterns package for rownames + set rownnames
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",package = "MutationalPatterns"))
df <- as.data.frame(df[,-1])
rownames(df) <- rownames(indel_counts)

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
    dplyr::select(all_of(samples)) %>%
    rowSums()
  
  plt <- plot_indel_contexts(df_types, condensed = TRUE, extra_labels = FALSE)
  plt <- plt + theme(axis.text.y = element_text(size=8), 
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 8),
                     axis.text.x = element_text(size = 6),
                     axis.title.y = element_blank(),
                     strip.background.y = element_blank(),
                     strip.text.y = element_blank(),
                     strip.background.x = element_blank(),
                     strip.text.x = element_blank(),
                     plot.margin = unit(c(0, 0, 0, 0), "in"))
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

figure <- plot_grid(plt1 + theme(legend.position = "none"), 
                    plt2 + theme(legend.position = "none"), 
                    plt3 + theme(legend.position = "none"),
                    plt4 + theme(legend.position = "none"),
                    plt5 + theme(legend.position = "none"),
                    plt6 + theme(legend.position = "none"),
                    plt7 + theme(legend.position = "none"), ncol = 1) 

pdf("results/figures/mutational_signatures/indel_mutation_context.pdf", width = 7.5, height = 9)
print(figure)
dev.off()

pdf("results/figures/mutational_signatures/indel_mutation_context_all.pdf", width = 7.5, height = 3)
print(plt1 + theme(legend.position = "none", strip.background.y = element_blank(), strip.text.y = element_blank()))
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
