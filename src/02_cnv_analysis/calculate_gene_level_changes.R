## ---------------------------
## Script Name: calculate_gene_level_changes.R
## Description: Calculate genes differentially impacted by CNV gains/losses
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-09-27
## ---------------------------

library(tidyverse)
library(optparse)

# PARSE COMMAND LINE OPTIONS
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "Path to data frame",
              metavar = "character"
  ),
  make_option(c("-m", "--metadata"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-g", "--gains"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-l", "--losses"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Load data
df <- readr::read_delim(opt$input, delim = "\t")
meta <- readxl::read_xlsx(opt$metadata, sheet = "tumors")

### Get number of metastatic and primary tumors
n_metastatic <- sum(meta$category == "Metastatic")
n_primary <- sum(meta$category == "Primary")

### Merge sample information onto cnv data
df <- df %>% dplyr::left_join(meta, by=c("TumorID" = "sample"))

### Select genes affected by gains and losses
gains <- df %>% dplyr::filter(Segment.Type %in% c("gain", "amp"))
losses <- df %>% dplyr::filter(Segment.Type %in% c("del", "loss"))

### Get count of genes affected by gains by tumor type
gains <- gains %>% dplyr::distinct(TumorID, Gene, category) %>% 
  dplyr::group_by(category, Gene) %>% dplyr::summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = category, values_from = count, values_fill = 0) %>%
  dplyr::mutate(Metastatic_wo = n_metastatic - Metastatic, Primary_wo = n_primary - Primary) %>%
  dplyr::mutate(MFreq = Metastatic / (Metastatic + Metastatic_wo), PFreq = Primary / (Primary + Primary_wo)) %>%
  dplyr::mutate(FreqDelta = abs(MFreq - PFreq))

### Calculate p.value and odds ratio
gains <- gains %>% dplyr::rowwise() %>%
  dplyr::mutate(p.value = fisher.test(as.table(rbind(c(Metastatic, Metastatic_wo), c(Primary, Primary_wo))))$p.value) %>%
  dplyr::mutate(odds.ratio = (Metastatic * Primary_wo) / (Primary * Metastatic_wo)) %>%
  dplyr::mutate(p.adjust = p.adjust(p.value, method="fdr"))

### Get count of genes affected by gains by tumor type
losses <- losses %>% dplyr::distinct(TumorID, Gene, category) %>% 
  dplyr::group_by(category, Gene) %>% dplyr::summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = category, values_from = count, values_fill = 0) %>%
  dplyr::mutate(Metastatic_wo = n_metastatic - Metastatic, Primary_wo = n_primary - Primary) %>%
  dplyr::mutate(MFreq = Metastatic / (Metastatic + Metastatic_wo), PFreq = Primary / (Primary + Primary_wo)) %>%
  dplyr::mutate(FreqDelta = abs(MFreq - PFreq))

### Calculate p.value and odds ratio
losses <- losses %>% dplyr::rowwise() %>%
  dplyr::mutate(p.value = fisher.test(as.table(rbind(c(Metastatic, Metastatic_wo), c(Primary, Primary_wo))))$p.value) %>%
  dplyr::mutate(odds.ratio = (Metastatic * Primary_wo) / (Primary * Metastatic_wo)) %>%
  dplyr::mutate(p.adjust = p.adjust(p.value, method="fdr"))

# Save results
readr::write_csv(losses, file = opt$losses)
readr::write_csv(gains, file = opt$gains)