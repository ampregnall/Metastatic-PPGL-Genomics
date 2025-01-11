## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-09
## ---------------------------


# Load Packages and Define Variables --------------------------------------

# Load packages
library(tidyverse)

# Define variables
read_depth <- 8
allele_freq <- 0.01
clin_var_filter <- c("Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic")

# Load Data ---------------------------------------------------------------

# List files
snv.files <- list.files("data/raw/snvs/", pattern = "*.csv", full.names = TRUE)

# Load into one data frame
snv.data <- map(snv.files, function(file) {
  read_delim(file, delim = ",",
             col_types = cols(
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
               VARSCAN2 = col_character()
             )
  )
}) %>% bind_rows()

n_mutations <- nrow(snv.data)
n_genes <- length(unique(snv.data$Gene))

# Save combined data
write_delim(snv.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.txt", delim = "\t")

# Filter by Allele Depth --------------------------------------------------

snv.data <- snv.data %>% filter(Tumor.AltDepth >= read_depth)
n_mutations_rd <- nrow(snv.data)
n_genes_rd <- length(unique(snv.data$Gene))

write_delim(snv.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.txt", delim = "\t")

# Filter by Clin.Var Annotations ------------------------------------------

snv.data <- snv.data %>% filter(!str_detect(tolower(ClinVar.SIG), "benign"))
n_mutations_clinvar <- nrow(snv.data)
n_genes_clinvar <- length(unique(snv.data$Gene))

# Apply filtering criteria to identify LOF/GOF mutations ------------------- 

snv.data <- snv.data %>% 
  filter(gnomAD.MAX_AF == "." | as.numeric(gnomAD.MAX_AF) < 0) %>%
  mutate(mutation.consequence = 
           case_when(str_detect(Variant.Consequence, "frameshift_variant|stop_gained") ~ "LOF",
                     Variant.Consequence %in% c("splice_acceptor_variant", "splice_donor_variant") ~ "LOF",
                     str_detect(Variant.Consequence, "missense_variant") & 
                       REVEL != "." & as.numeric(REVEL) > 0.5 ~ "LOF",
                     str_detect(Variant.Consequence, "missense_variant") & 
                       (ClinVar.SIG %in% clin_var_filter | grepl("COSV", Existing.variation)) ~ "GOF",
                     TRUE ~ "Unknown"))

lof.data <- snv.data %>% filter(mutation.consequence == "LOF")
gof.data <- snv.data %>% filter(mutation.consequence == "GOF")

n_mutations_lof <- nrow(lof.data)
n_genes_lof <- length(unique(lof.data$Gene))
n_mutations_gof <- nrow(gof.data)
n_genes_gof <- length(unique(gof.data$Gene))

# Save filtered data
write_delim(lof.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.txt", delim = "\t")
write_delim(gof.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.txt", delim = "\t")
         
# Filter by COSMIC --------------------------------------------------------

cgc <- read_delim("metadata/cancer_census_genes_all_v98.tsv", delim = "\t")

lof.data <- lof.data %>% filter(Gene %in% cgc$`Gene Symbol`)
gof.data <- gof.data %>% filter(Gene %in% cgc$`Gene Symbol`) 

n_mutations_lof_cgc <- nrow(lof.data)
n_genes_lof_cgc <- length(unique(lof.data$Gene))
n_mutations_gof_cgc <- nrow(gof.data)
n_genes_gof_cgc <- length(unique(gof.data$Gene))

write_delim(lof.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.LOF.CGC.txt", delim = "\t")
write_delim(gof.data, "data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.GOF.CGC.txt", delim = "\t")
