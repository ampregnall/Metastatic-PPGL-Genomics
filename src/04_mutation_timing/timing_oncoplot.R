## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-21
## Copyright (c) Andrew Pregnall, 2025
## ---------------------------

# Load packages
#library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Define variables
arms <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
          "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q",
          "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q",
          "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q")

col <- c("#377EB8",
         "#4DAF4A",
         "#984EA3", 
         "#E41A1C",
         "#999999")



# Load required data ------------------------------------------------------

# Load cytoband coordinates
cytos <- readr::read_delim("metadata/GrCh38.arm.coordinates.txt")

# Load sample information
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
meta <- meta %>% dplyr::arrange(order)
samples <- meta %>% dplyr::select(sample)

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
}) %>% bind_rows() 

# Get missing copy number information from SNV files
cnv.info <- snvs %>% 
  distinct(Tumor.ID, CNID, .keep_all = TRUE) %>%
  dplyr::select(Tumor.ID, CNID, chromosome, start.pos, end.pos)

# Combine data -----
cnvs <- cnvs %>% left_join(cnv.info, by = c("Tumor.ID", "CNID"))
cnvs <- cnvs %>% left_join(meta[,c("sample", "wholeGenomeDuplication")], by = c("Tumor.ID" = "sample"))

# Identify chromosomal arm of copy number variants ------------------------

cnvs <- cnvs %>% 
  left_join(cytos, by = c("chromosome" = "chrom")) 

cnvs <- cnvs %>%
  # Create rows for both conditions at once with case_when
  filter(start.pos >= p_start & start.pos <= p_end & end.pos > q_start & end.pos <= q_end) %>%
  mutate(
    start.pos = case_when(
      end.pos > q_start & end.pos <= q_end ~ q_start,
      TRUE ~ start.pos
    ),
    end.pos = case_when(
      start.pos >= p_start & start.pos <= p_end ~ p_end,
      TRUE ~ end.pos
    )
  ) %>%
  # Combine with the filtered original dataset
  dplyr::bind_rows(
    cnvs %>%
      filter((start.pos >= p_start & end.pos <= p_end) | (start.pos >= q_start & end.pos <= q_end))
  )


cnvs <- cnvs %>% 
  mutate(arm = case_when(start.pos >= p_start & end.pos <= p_end ~ str_c(chromosome, "p"),
                         start.pos >= q_start & end.pos <= q_end ~ str_c(chromosome, "q"))) %>% 
  mutate(arm = factor(arm, levels = arms))
  
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
      dplyr::select(arm, time, time.up, time.lo, n.snv_mnv) %>%
      group_by(arm) %>%
      aggregate_segments() %>%
      filter(arm %in% arms) %>%
      arrange(arm)
  } else {
    # Whole-genome duplication scenario
    chrom_result <- bb %>%
      as_tibble() %>%
      dplyr::select(arm, time, time.up, time.lo, n.snv_mnv) %>%
      group_by(arm) %>%
      aggregate_segments() %>%
      filter(arm %in% arms) %>%
      arrange(arm)
    
    wgd_result <- bb %>%
      dplyr::select(time, time.up, time.lo, n.snv_mnv) %>%
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
  })

# Clonal/subclonal information --------------------------------------------

# Count number of clonal/subclonal mutations by sample
clonality.summary <- snvs %>% 
  group_by(Tumor.ID, CLS) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = CLS, values_from = count, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ . / rowSums(across(where(is.numeric)))))

clonality.summary <- as.data.frame(clonality.summary)
rownames(clonality.summary) = clonality.summary[, 1]
clonality.summary <- clonality.summary[match(samples$sample, rownames(clonality.summary)), ]
clonality.summary = clonality.summary[, -1]
clonality.summary <- clonality.summary[rownames(clonality.summary) %in% meta$sample,]

# Aggregate timing information by chromosomal arm -------------------------

timing <- cnvs.agg %>%
  dplyr::select(sample, arm, time) %>%
  pivot_wider(names_from = arm, values_from = time) 

timing <- left_join(samples, timing, by = c("sample" = "sample")) 
timing <- as.data.frame(timing)
rownames(timing) = timing[, 1]
timing <- timing[match(samples$sample, rownames(timing)), ]
timing = timing[, -1]

### Add missing chromosomes
timing$`13p` <- NA
timing$`14p` <- NA
timing$`15p` <- NA
timing$`22p` <- NA
timing <- timing %>% dplyr::rename(WGD = "NA")


### Define colors for Oncoplot
col_fun = colorRamp2(c(0, 0.5, 1), c("#3c7b4a", "white", "#6b3e8a"))

col_category <- c("Primary" = "#2D3A7F", 
                  "Metastatic" = "lightgrey")

col_tumor <- c("PCC" = "#879eb3", 
               "PGL"="#542E71")

col_germline <- c("SDHA" = "#D3672D", 
                  "SDHB" = "#1B4367", 
                  "SDHC" = "#5F8528", 
                  "RET" =  "#A9361E", 
                  "WT" =   "#3F8B87")


# Create plot
plot <- Heatmap(as.matrix(timing),
  col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE, column_order = c(arms, "WGD"),
  gap = 1, show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = TRUE,
  left_annotation = rowAnnotation(Tumor = anno_simple(x = meta$tumor_type, col = col_tumor),
                                  Type = anno_simple(x = meta$category, col = col_category),
                                  Germline = anno_simple(meta$germline, col = col_germline), show_legend = FALSE),
  right_annotation = rowAnnotation(Clonality = anno_barplot(clonality.summary,
    gp = gpar(fill = col, border = FALSE), width = unit(1, "in"), border = FALSE
  ))
)

### Save results
pdf("results/figures/mutational_timing/timing_oncoplot.pdf", width = 15, height = 7.5)
print(plot)
dev.off()