library(signature.tools.lib)
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
  make_option(c("-o", "--output"),
              type = "character",
              help = "output file name [default = %default]",
              metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

df <- readr::read_delim(opt$input, delim = "\t")

df_snv <- df %>% dplyr::filter(Variant.Class == "SNV") %>%
  dplyr::select(Tumor.ID, Chr, Start, REF, ALT) %>%
  dplyr::rename(chr = Chr, position = Start)

df_snv_split <- df_snv %>% base::split(f = as.factor(.$Tumor.ID))

df_snv_split <- lapply(df_snv_split, signature.tools.lib::tabToSNVcatalogue, genome.v = "hg38")

df_snv_catalogs <- purrr::map(df_snv_split, 2) %>%
  purrr::map(1) %>% as_tibble()

# SAVE RESULTS
readr::write_delim(df_snv_catalogs, opt$output, delim = ",")