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
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Read in GOF candidate variants
df <- readr::read_delim(opt$input)

### Filter based on ClinVar or COSMIC annotation 
clin_var_filter <- c("Uncertain_significance", "Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic")
df1 <- df %>% dplyr::filter(ClinVar.SIG %in% clin_var_filter | grepl("COSV", Existing.variation))

### Save results
readr::write_delim(df1, opt$output, delim = "\t")