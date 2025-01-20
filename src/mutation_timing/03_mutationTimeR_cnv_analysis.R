## ---------------------------
## Script Name: 03_mutationTimeR_cnv_analysis.R
## Description: Analyze CNV data from MutationTimeR
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-14
## ---------------------------

# Load packages
library(tidyverse)
library(readxl)

# Define variables
prgn <- RColorBrewer::brewer.pal(11,"PRGn")
set1 <- RColorBrewer::brewer.pal(9,"Set1")
col <- colorRampPalette(set1[c(4,9,3)])(10)
chr <- c(seq(1:22), "X")

# Get WGD samples
meta <- read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors") 

# Load data ---------------------------------------------------------------

# Load CNV data with timing estimates
cnv.files <- list.files(
  path = "data/processed/mutation_timing/MutationTimeR/timing",
  pattern = "*.csv", full.names = FALSE, recursive = FALSE
)

cnvs <- map(cnv.files, function(file) {
  tmp <- read_csv(paste0("data/processed/mutation_timing/MutationTimeR/timing/", file),
    col_names = c("CNID", "type", "time", "time.lo", "time.up", "time.2nd", 
                  "time.2nd.lo", "time.2nd.up", "time.star", "n.snv_mnv"),
    col_types = cols(
      CNID = col_double(),
      type = col_character(),
      time = col_double(),
      time.lo = col_double(),
      time.up = col_double(),
      time.2nd = col_double(),
      time.2nd.lo = col_double(),
      time.2nd.up = col_double(),
      time.star = col_character(),
      n.snv_mnv = col_double()
    ), skip = 1
  )
  
  # Add sample name
  tmp <- tmp %>% mutate(Tumor.ID = str_remove(file, "_mT.csv"))
  return(tmp)
}) %>% bind_rows()

# Remove rows with NA timing estimates
cnvs <- cnvs %>% filter(!is.na(time))

# Load variant information that has CNID linked to chromosome
files <- list.files(
  path = "data/processed/mutation_timing/MutationTimeR/variant_timing",
  pattern = "*.csv", full.names = FALSE, recursive = FALSE
)

snvs <- map(files, function(file) {
  read_csv(paste0("data/processed/mutation_timing/MutationTimeR/variant_timing/", file),
           col_types = cols(
             ...1 = col_double(),
             Tumor.ID = col_character(),
             Normal.ID = col_character(),
             Chr = col_character(),
             Start = col_double(),
             REF = col_character(),
             ALT = col_character(),
             FILTER = col_character(),
             Gene = col_character(),
             Gene.Accession = col_character(),
             Variant.Class = col_character(),
             Variant.Type = col_character(),
             Variant.Consequence = col_character(),
             HGVSc = col_character(),
             HGVSp = col_character(),
             Feature.Type = col_character(),
             Feature.Accession = col_character(),
             Bio.type = col_character(),
             Existing.variation = col_character(),
             EXON = col_character(),
             INTRON = col_character(),
             STRAND = col_double(),
             cDNA.position = col_character(),
             CDS.position = col_character(),
             Protein.position = col_character(),
             Amino.acids = col_character(),
             Codons = col_character(),
             SpliceAI.DS_AG = col_character(),
             SpliceAI.DS_AL = col_character(),
             SpliceAI.DS_DG = col_character(),
             SpliceAI.DS_DL = col_character(),
             REVEL = col_character(),
             gnomAD.AF = col_character(),
             gnomAD.AFR = col_character(),
             gnomAD.AMR = col_character(),
             gnomAD.ASJ = col_character(),
             gnomAD.EAS = col_character(),
             gnomAD.FIN = col_character(),
             gnomAD.NFE = col_character(),
             gnomAD.OTH = col_character(),
             gnomAD.SAS = col_character(),
             gnomAD.MAX_AF = col_character(),
             gnomAD.MAX_POPS = col_character(),
             LOFTEE.lof = col_character(),
             LOFTEE.filter = col_character(),
             LOFTEE.flags = col_character(),
             LOFTEE.info = col_character(),
             ClinVar = col_character(),
             ClinVar.SIG = col_character(),
             ClinVar.REVSTAT = col_character(),
             ClinVar.DN = col_character(),
             Tumor.Zyg = col_character(),
             Tumor.Depth = col_double(),
             Tumor.AltDepth = col_double(),
             Tumor.AltFrac = col_double(),
             Normal.Zyg = col_character(),
             Normal.Depth = col_double(),
             Normal.AltDepth = col_double(),
             Normal.AltFrac = col_double(),
             LANCET = col_character(),
             MUTECT2 = col_character(),
             STRELKA2 = col_character(),
             VARDICT = col_character(),
             VARSCAN2 = col_character(),
             MutCN = col_double(),
             MutDeltaCN = col_double(),
             MajCN = col_double(),
             MinCN = col_double(),
             MajDerCN = col_double(),
             MinDerCN = col_double(),
             CNF = col_double(),
             CNID = col_double(),
             pMutCN = col_double(),
             pGain = col_double(),
             pSingle = col_double(),
             pSub = col_double(),
             pMutCNTail = col_double(),
             CLS = col_character(),
             CN = col_double(),
             chromosome = col_character(),
             start.pos = col_double(),
             end.pos = col_double(),
             Bf = col_double(),
             N.BAF = col_double(),
             sd.BAF = col_double(),
             depth.ratio = col_double(),
             N.ratio = col_double(),
             sd.ratio = col_double(),
             CNt = col_double(),
             A = col_double(),
             B = col_double(),
             LPP = col_double()
           )
  )
}) %>% bind_rows() %>% 
  distinct(Tumor.ID, CNID, .keep_all = TRUE) %>%
  select(Tumor.ID, CNID, chromosome, start.pos, end.pos)

# Combine data -----
cnvs <- cnvs %>% left_join(snvs, by = c("Tumor.ID", "CNID"))
cnvs <- cnvs %>% left_join(meta[,c("sample", "wholeGenomeDuplication")], by = c("Tumor.ID" = "sample"))

# Aggregate timing data per chromosome -----

aggregatePerChromosome <- function(bb, isWgd = FALSE) {
  
  # Helper function for aggregation
  aggregate_segments <- function(df) {
    df %>%
      summarise(
        time = weighted.mean(time, n.snv_mnv, na.rm = TRUE),
        n = sum(n.snv_mnv[!is.na(time)], na.rm = TRUE),
        sd = sd(time, na.rm = TRUE),
        ci = weighted.mean(time.up - time.lo, n.snv_mnv, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Main aggregation logic
  if (!isWgd) {
    # Group by chromosome and aggregate
    result <- bb %>%
      select(chromosome, time, time.up, time.lo, n.snv_mnv) %>%
      group_by(chromosome) %>%
      aggregate_segments() %>%
      filter(chromosome %in% c(as.character(1:22), "X")) %>%
      arrange(chromosome)
  } else {
    # Whole-genome duplication scenario
    chrom_result <- bb %>%
      as_tibble() %>%
      select(chromosome, time, time.up, time.lo, n.snv_mnv) %>%
      group_by(chromosome) %>%
      aggregate_segments() %>%
      filter(chromosome %in% c(as.character(1:22), "X")) %>%
      arrange(chromosome)
    
    wgd_result <- bb %>%
      select(time, time.up, time.lo, n.snv_mnv) %>%
      aggregate_segments()
    
    result <- chrom_result %>%
      bind_rows(wgd_result %>% mutate(seqnames = "WGD"))
  }
  
  return(result)
}

cnvs.agg <- cnvs %>%
  group_by(Tumor.ID, wholeGenomeDuplication) %>%
  group_split() %>%
  map_dfr(~ {
    aggregated <- aggregatePerChromosome(.x, isWgd = unique(.x$wholeGenomeDuplication))
    aggregated %>%
      mutate(sample = unique(.x$Tumor.ID))
  }) %>% mutate(chromosome = factor(chromosome, levels = chr))


# Plot timing estimates across all samples ------------------------------

cnvs.agg.nWGD <- cnvs.agg %>% filter(is.na(seqnames))
cnvs.agg.WGD <- cnvs.agg %>% filter(!is.na(seqnames))

# Plot distribution of timing estimates aggregated across chromosomes and samples
plt1 <- ggplot(cnvs.agg.nWGD, aes(x = time, fill = ..x..)) +
  geom_histogram(binwidth = 0.1, position = "identity") +
  scale_fill_gradientn(colors = rev(prgn), limits = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Timing Estimate", y = "Count", fill = "Sample ID") +
  ggtitle("Distribution of Timing Estimates Across All Samples") +
  theme(legend.position = "none")

pdf("results/figures/mutational_timing/mutation_timing_cnv_agg.pdf", width = 10, height = 5)
print(plt1)
dev.off()

# Plot distribution of WGD timing estimates aggregated samples
plt2 <- ggplot(cnvs.agg.WGD, aes(x = time, fill = ..x..)) +
  geom_histogram(binwidth = 0.1, position = "identity") +
  scale_fill_gradientn(colors = rev(prgn), limits = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Timing Estimate", y = "Count", fill = "Sample ID") +
  ggtitle("Distribution of WGD Timing Estimates Across All Samples") +
  theme(legend.position = "none")

pdf("results/figures/mutational_timing/mutation_timing_wgd_agg.pdf", width = 10, height = 5)
print(plt2)
dev.off()

# Plot distribution of timing estimates aggregated across samples, split by chromosome
plt3 <- ggplot(cnvs.agg.nWGD, aes(x = time, fill = ..x..)) +
  geom_histogram(binwidth = 0.1, position = "identity") +
  scale_fill_gradientn(colors = rev(prgn), limits = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(x = "Timing Estimate", y = "Count", fill = "Sample ID") +
  ggtitle("Distribution of Timing Estimates Across All Samples") +
  theme(legend.position = "none") +
  facet_wrap(~chromosome, scales = "free_y")

pdf("results/figures/mutational_timing/mutation_timing_cnvs_by_chr.pdf", width = 10, height = 5)
print(plt3)
dev.off()