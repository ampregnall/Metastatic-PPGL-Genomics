library(tidyverse)
library(optparse)
source("src/mutation_timing/estimate_mutation_timing.R")

option_list <- list(
  make_option(c("-s", "--snv"),
              type = "character",
              default = NULL,
              help = "Path to data frame containing somatic variants"),
  make_option(c("-c", "--cnv"),
              type = "character",
              default = NULL,
              help = "Path to data frame containing copy number variants"),
  make_option(c("-m", "--metadata"),
              type = "character",
              default = NULL,
              help = "Path to data frame containing sample purity estimates"),
  make_option(c("--sample"), 
              type = "character",
              default = NULL,
              help = "Sample name"),
  make_option(c("-o", "--output"),
              type = "character",
              default = "out.txt",
              help = "Output file name: [default = %default]"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Load somatic and copy number variants 
copy_number_vars <- readr::read_csv(opt$cnv)
somatic_vars <- readr::read_csv(opt$snv)

### Load purity and ploidy estimates for samples
estimates <- readxl::read_xlsx(opt$metadata) 
purity <- estimates[estimates$sample == opt$sample,]$cellularity

### Create unique mutation ID and select distinct mutations
mutations <- somatic_vars %>% dplyr::mutate(MutID = stringr::str_c(Chr, Start, Gene, REF, ALT, sep = "-")) %>%
  dplyr::select(Tumor.ID, Chr, Start, Gene, REF, ALT, MutID, Tumor.Depth, Tumor.AltDepth, Tumor.AltFrac, Gene.Accession) %>%
  dplyr::distinct(MutID, .keep_all = TRUE) %>%
  dplyr::mutate(Chr = as.character(Chr))

### Extract copy number info for each mutation; remove variants on X chromosome
cnv_per_mutation <- mutations %>% 
  dplyr::left_join(copy_number_vars, by=c("Gene"="Gene...2", "Chr"="SV.Chrom"), relationship = "many-to-many") %>%
  dplyr::filter(Start >= SV.Start & Start <= SV.End) %>% dplyr::filter(Chr != "X")

### Extract max copy number
max_cn <- max(cnv_per_mutation$Segment.CN)

### Estimate 95% CI of VAF
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(Tumor.AltFrac.95 = prop.test(Tumor.AltDepth, Tumor.Depth)$conf[2]) %>%
  dplyr::mutate(Tumor.AltFrac.05 = prop.test(Tumor.AltDepth, Tumor.Depth)$conf[1])

### Estimate mutation multiplicity in samples
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(Mutation.Multiplicity = (Tumor.AltFrac / purity) * ((purity * Segment.CN)+ 2 * (1 - purity))) %>%
  dplyr::mutate(Mutation.Multiplicity.95 = (Tumor.AltFrac.95 / purity) * ((purity * Segment.CN)+ 2 * (1 - purity))) %>%
  dplyr::mutate(Mutation.Multiplicity.05 = (Tumor.AltFrac.05 / purity) * ((purity * Segment.CN) + 2 * (1 - purity)))

estimate_cancer_cell_fraction <- function(n_tumor_alt, n_tumor, CNt) {
  # TODO: move function into own file? 
  # TODO: add more documentation below to explain steps
  # Calculate distribution of VAFs assuming CCF of x: (x*purity) / (2*(1-purity)+purity*CNt
  # Calculate probability of each CCF given tumor/total reads using binomial probability distribution
  if (CNt == 0) {
    results <- list(NA, NA, NA, NA, NA)
    names(results) <- c("CCF", "CCF.05", "CCF.95", "Prob.Clonal", "Prob.Subclonal")
    return(results)
  } 
  else {
    x <- dbinom(n_tumor_alt, n_tumor, 
                prob = sapply(seq(0.01, 1, 0.01), function(x) (x*purity) / (2*(1-purity)+purity*CNt)))
    
    if(min(x) == 0) {
      x[length(x)] <- 1}
    names(x) <- seq(0.01, 1, 0.01)
    
    # Normalize probability distribution
    x_norm <- x / sum(x) 
    
    
    x_sort <- sort(x_norm, decreasing = TRUE)
    x_cumalative_likelihood <- cumsum(x_sort)
    
    n = sum(x_cumalative_likelihood < 0.95) + 1
    threshold <- x_sort[n]
    conf_int <- x[x_norm >= threshold]
    cancer_cell_fractions <- as.numeric(names(conf_int))
    ccf <- cancer_cell_fractions[which.max(conf_int)]
    ccf.05 <- cancer_cell_fractions[1]
    ccf.95 <- cancer_cell_fractions[length(cancer_cell_fractions)]
    prob_subclonal <- sum(x_norm[1:90])
    prob_clonal <- sum(x_norm[91:100])
    results <- list(ccf, ccf.05, ccf.95, prob_clonal, prob_subclonal)
    names(results) <- c("CCF", "CCF.05", "CCF.95", "Prob.Clonal", "Prob.Subclonal")
    return(results)
  }
}

### Calculate CCF for the sample
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(results = list(estimate_cancer_cell_fraction(Tumor.AltDepth, Tumor.Depth, Segment.CN))) %>%
  tidyr::unnest_wider(results) %>%
  dplyr::ungroup()

### Estimate mutation timing
cnv_per_mutation <- cnv_per_mutation %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(timing = estimate_mutation_timing(max_cn = max_cn, 
                                                  major_cn = Segment.CNa, 
                                                  n_tumor = Tumor.Depth, 
                                                  vaf = Tumor.AltFrac, 
                                                  purity = purity))

### SAVE RESULTS
readr::write_delim(cnv_per_mutation, opt$output, delim="\t")