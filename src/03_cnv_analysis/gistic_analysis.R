## ---------------------------
## Script Name: 01_gistic_analysis.R
## Description: Create GISTIC and IGV input. Plot results of GISTIC analysis 
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-24
## Copyright (c) Andrew Pregnall, 2025
## ---------------------------

# Load packages
library(purrr)
library(stringr)
library(dplyr)
library(readr)
library(maftools)

# Define variables
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
prim <- filter(meta, category == "Primary")$sample
met <- filter(meta, category == "Metastatic")$sample

# Load data -----------------------------------------------------------

# Load all files in directory
files <- list.files(
  path = "data/raw/cnvs/sequenza", pattern = "*.txt", full.names = FALSE, recursive = FALSE
)

# Read all CSV files into a list
data_list <- map(files, function(file) {
  result <- read_delim(paste0("data/raw/cnvs/sequenza/", file) , delim = "\t")
  result %>% mutate(sample = stringr::str_split(file, "_")[[1]][1])
})

# Combine all data frames
df <- bind_rows(data_list)

# Preprocess data for GISTIC ----------------------------------------------

# Prepare data for GISTIC analysis
df <- df %>% 
  select(sample, chromosome, start.pos, end.pos, N.BAF, depth.ratio) %>%
  distinct(sample, chromosome, start.pos, end.pos, N.BAF, depth.ratio) %>% 
  mutate(Seg.CN = log2(2 * depth.ratio) - 1) %>%
  select(sample, chromosome, start.pos, end.pos, N.BAF, Seg.CN)

# Subset data
df.prim <- df %>% filter(sample %in% prim)
df.met <- df %>% filter(sample %in% met)

# Save data
write_delim(df, "data/processed/cnvs/gistic/mPPGL.gistic.input.all.txt", delim = "\t")
write_delim(df, "data/processed/cnvs/gistic/mPPGL.gistic.input.all.seg", delim = "\t")
write_delim(df.prim, "data/processed/cnvs/gistic/mPPGL.gistic.input.primaries.txt", delim = "\t")
write_delim(df.met, "data/processed/cnvs/gistic/mPPGL.gistic.input.metastatic.txt", delim = "\t")


# RUN GISTIC ------------------------------------------------------------------

## ---------------------------
## Output files above were used to run GISTIC using the GenePattern server
## Results were manually saved to output directors below
## Seg file above was loaded into IGV to create visualization of CNVs
## ---------------------------

# GISTIC plotting -------------------------------------------------------

# Load GISTIC results
gistic.met = readGistic(gisticAllLesionsFile = "data/processed/cnvs/gistic/gistic_metastatic/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "data/processed/cnvs/gistic/gistic_metastatic/amp_genes.conf_90.txt",
                         gisticDelGenesFile = "data/processed/cnvs/gistic/gistic_metastatic/del_genes.conf_90.txt",
                         gisticScoresFile = "data/processed/cnvs/gistic/gistic_metastatic/scores.gistic")

gistic.prim = readGistic(gisticAllLesionsFile = "data/processed/cnvs/gistic/gistic_primary/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "data/processed/cnvs/gistic/gistic_primary/amp_genes.conf_90.txt",
                         gisticDelGenesFile = "data/processed/cnvs/gistic/gistic_primary/del_genes.conf_90.txt",
                         gisticScoresFile = "data/processed/cnvs/gistic/gistic_primary/scores.gistic")

# Save results
pdf("results/figures/cnv_analysis/gistic_coplot_amp.pdf", width = 9, height = 12)
coGisticChromPlot(gistic1 = gistic.prim, 
                  gistic2 = gistic.met, 
                  g1Name = "Primary", 
                  g2Name = "Metastatic", 
                  type = 'Amp')
dev.off()

pdf("results/figures/cnv_analysis/gistic_coplot_del.pdf", width = 9, height = 12)
coGisticChromPlot(gistic1 = gistic.prim, 
                  gistic2 = gistic.met, 
                  g1Name = "Primary", 
                  g2Name = "Metastatic", 
                  type = 'Del', )
dev.off()

# Load CGC data and filter GISTIC results for annotations
cgc <- read_delim("metadata/cancer_census_genes_all_v98.tsv", delim = "\t")
cgc <- filter(cgc, Tier == 1)
met.genes <- gistic.met@gene.summary %>% filter(Hugo_Symbol %in% cgc$'Gene Symbol')
prim.genes <- gistic.prim@gene.summary %>% filter(Hugo_Symbol %in% cgc$'Gene Symbol')
