library(purrr)
library(stringr)
library(dplyr)
library(readr)

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

# Prepare data for GISTIC analysis
df <- df %>% dplyr::select(sample, chromosome, start.pos, end.pos, N.BAF, depth.ratio) %>%
  dplyr::distinct(sample, chromosome, start.pos, end.pos, N.BAF, depth.ratio) %>% 
  dplyr::mutate(Seg.CN = log2(2 * depth.ratio) - 1) %>%
  dplyr::select(sample, chromosome, start.pos, end.pos, N.BAF, Seg.CN)

# Load metadata
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
prim <- dplyr::filter(meta, category == "Primary")$sample
met <- dplyr::filter(meta, category == "Metastatic")$sample

# Subset data
df.prim <- df %>% filter(sample %in% prim)
df.met <- df %>% filter(sample %in% met)

# Save data
write_delim(df, "data/processed/cnvs/gistic/mPPGL.gistic.input.all.txt", delim = "\t")
write_delim(df.prim, "data/processed/cnvs/gistic/mPPGL.gistic.input.primaries.txt", delim = "\t")
write_delim(df.met, "data/processed/cnvs/gistic/mPPGL.gistic.input.metastatic.txt", delim = "\t")

# GISTIC plotting -------------------------------------------------------

library(maftools)
source("src/utilities/gistic_coplot.R")
source("src/utilities/readSegs.R")

# Load GISTIC results
gistic.all = readGistic(gisticAllLesionsFile = "data/processed/cnvs/gistic/gistic_all/all_lesions.conf_90.txt", 
                        gisticAmpGenesFile = "data/processed/cnvs/gistic/gistic_all/amp_genes.conf_90.txt",
                        gisticDelGenesFile = "data/processed/cnvs/gistic/gistic_all/del_genes.conf_90.txt",
                        gisticScoresFile = "data/processed/cnvs/gistic/gistic_all/scores.gistic")

gistic.met = readGistic(gisticAllLesionsFile = "data/processed/cnvs/gistic/gistic_metastatic/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "data/processed/cnvs/gistic/gistic_metastatic/amp_genes.conf_90.txt",
                         gisticDelGenesFile = "data/processed/cnvs/gistic/gistic_metastatic/del_genes.conf_90.txt",
                         gisticScoresFile = "data/processed/cnvs/gistic/gistic_metastatic/scores.gistic")

gistic.prim = readGistic(gisticAllLesionsFile = "data/processed/cnvs/gistic/gistic_primary/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile = "data/processed/cnvs/gistic/gistic_primary/amp_genes.conf_90.txt",
                         gisticDelGenesFile = "data/processed/cnvs/gistic/gistic_primary/del_genes.conf_90.txt",
                         gisticScoresFile = "data/processed/cnvs/gistic/gistic_primary/scores.gistic")

# Create overall plot
gisticChromPlot(gistic = gistic.all, markBands = "all")

coGisticChromPlot(gistic1 = gistic.prim, gistic2 = gistic.met, g1Name = "Primary", g2Name = "Metastatic", type = 'Amp')
coGisticChromPlot(gistic1 = gistic.prim, gistic2 = gistic.met, g1Name = "Primary", g2Name = "Metastatic", type = 'Del')



