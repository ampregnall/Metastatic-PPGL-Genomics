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
  make_option(c("-c", "--cosmic"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-o", "--output"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  )
)

### Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Load input data
df <- readr::read_delim(opt$input)
cgc <- readr::read_delim(opt$cosmic)

### Select for COSMIC genes in GOF candidates
df1 <- df %>% dplyr::filter(Gene %in% cgc$`Gene Symbol`)
df1 <- df1 %>% dplyr::left_join(dplyr::select(cgc, `Gene Symbol`, `Role in Cancer`), by=c("Gene"="Gene Symbol"))

### Save results
readr::write_delim(df1, opt$output, delim = "\t")
