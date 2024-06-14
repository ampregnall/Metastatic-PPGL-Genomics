library(tidyverse)

### Load somatic and copy number variants 
copy_number_vars <- readr::read_delim("data/raw/cnvs/sequenza/PP002-DZ3A_segments.txt", delim = "\t")
somatic_vars <- readr::read_csv("data/raw/snvs/PP002-DZ3A.ffpe_scenario.local-fdr.paired_somatic.vep.report.csv")

### Load purity and ploidy estimates for samples
# estimates <- readxl::read_xlsx(opt$metadata) 
estimates <- readxl::read_xlsx("metadata/Sequenza-Ploidy-Estimates.xlsx") 
purity <- estimates[estimates$sample == "PP002-DZ3A",]$cellularity
ploidy <- estimates[estimates$sample == "PP002-DZ3A",]$ploidy

### Create unique mutation ID and select distinct mutations
mutations <- somatic_vars %>% dplyr::mutate(MutID = stringr::str_c(Chr, Start, Gene, REF, ALT, sep = "-")) %>%
  dplyr::select(Chr, Start, Gene, REF, ALT, MutID, Tumor.Depth, Tumor.AltDepth, Tumor.AltFrac, Gene.Accession) %>%
  dplyr::distinct(MutID, .keep_all = TRUE) %>%
  dplyr::mutate(Chr = as.character(Chr))

### Extract copy number info for each mutation
cnv_per_mutation <- mutations %>% dplyr::left_join(copy_number_vars, by=c("Chr"="chromosome"), relationship = "many-to-many") %>%
  dplyr::filter(Start >= start.pos & Start <= end.pos) 

### Estimate 95% CI of VAF
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(Tumor.AltFrac.95 = prop.test(Tumor.AltDepth, Tumor.Depth)$conf[2]) %>%
  dplyr::mutate(Tumor.AltFrac.05 = prop.test(Tumor.AltDepth, Tumor.Depth)$conf[1])

### Estimate mutation multiplicity in samples
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(Mutation.Multiplicity = (Tumor.AltFrac / purity) * ((purity * CNt)+ 2 * (1 - purity))) %>%
  dplyr::mutate(Mutation.Multiplicity.95 = (Tumor.AltFrac.95 / purity) * ((purity * CNt)+ 2 * (1 - purity))) %>%
  dplyr::mutate(Mutation.Multiplicity.05 = (Tumor.AltFrac.05 / purity) * ((purity * CNt) + 2 * (1 - purity)))

### Get maximum copy number from sample
max_cn <- max(cnv_per_mutation$CNt)

estimate_cancer_cell_fraction <- function(n_tumor_alt, n_tumor, CNt) {
  
  # Calculate distribution of VAFs assuming CCF of x: (x*purity) / (2*(1-purity)+purity*CNt
  # Calculate probability of each CCF given tumor/total reads using binomial probability distribution
  x <- dbinom(n_tumor_alt, n_tumor, prob = sapply(seq(0.01, 1, 0.01), function(x) (x*purity) / (2*(1-purity)+purity*CNt)))
  
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
  results
}

### Calculate CCF for the sample
cnv_per_mutation <- cnv_per_mutation %>% dplyr::rowwise() %>%
  dplyr::mutate(results = list(estimate_cancer_cell_fraction(Tumor.AltDepth, Tumor.Depth, CNt))) %>%
  tidyr::unnest_wider(results) %>%
  dplyr::ungroup()

# SANDBOX -----------------------------------------------------------------

### Define function to create combinations of A/B alleles; borrowed from Sequenza 3.0
mufreq.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
  cn_ratio_vect <- seq(from = CNt.min / CNn, to = CNt.max / CNn, by = 1 / CNn)
  CNt <- cn_ratio_vect * CNn
  mut_comb <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
  times_muts <- sapply(mut_comb, length)
  data.frame(CNn = CNn, CNt = rep(CNt, times = times_muts),
             Mt = unlist(mut_comb))
}

### Define matrix; filter for at least 1 B allele; calculate theoretical frequencies
types <- mufreq.types.matrix(CNt.min = 1, CNt.max = max_cn, CNn = 2)
types <- types %>% dplyr::filter(Mt >= 1)

### Define function to calculate theoretical mutation frequency at various copy number states
theoretical.mufreq <- function(Mt, CNt, CNn = 2, cellularity) {
  normal.alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
  all.alleles    <- (CNt * cellularity) + CNn * (1 - cellularity)
  1 - (normal.alleles / all.alleles)
}

### Calculate expected mutation frequency at various copy number states
types <- types %>% dplyr::rowwise() %>%
  dplyr::mutate(test = theoretical.mufreq(Mt = Mt, CNt, CNn, cellularity = purity))


mufreq.dpois <- function(mufreq, mufreq.model, depth.t, seq.errors = 0.01) {
  mufreq.model[mufreq.model == 0] <- seq.errors
  n.success <- round(mufreq * depth.t, 0)
  dpois(x = n.success, lambda = mufreq.model * depth.t)
}

get.Mt <- function(F, depth.t, types, CNt, CNn, Mt){
  types <- types[types$CNn == CNn, ]
  l <- mufreq.dpois(mufreq = F, types$F[types$CNt== CNt&types$Mt<=Mt], depth.t = depth.t)
  l <- l/sum(l)
  L <- data.frame(l = l, Mt = types$Mt[types$CNt== CNt&types$Mt<=Mt])
  return(L)
}

