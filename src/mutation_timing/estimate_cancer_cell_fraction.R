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

