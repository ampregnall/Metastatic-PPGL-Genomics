## ---------------------------
## Script Name: jaccard_similarity_index.R
## Description: Calculate JSI based on Hu et al (Nat Gen 2020)
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-14
## ---------------------------t

# Load packages
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)


# Define variables
clonal <- c("clonal [late]", "clonal [NA]", "clonal [early]")
subclonal <- c("subclonal")

# Load data -------
files <- list.files(
  path = "data/processed/mutation_timing/MutationTimeR/variant_timing",
  pattern = "*.csv", full.names = FALSE, recursive = FALSE
)

df <- map(files, function(file) {
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

# Add unique mutation ID
df <- df %>% mutate(MutID = str_c(Chr, Start, Gene, REF, ALT, sep = "-"))

# Load tumor pairings
jsi <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "paired")

# Define empty vectors to store mutation information
metastasis_private_clonal <- vector(mode = "numeric", length = nrow(jsi))
primary_private_clonal <- vector(mode = "numeric", length = nrow(jsi))
metastasis_private_subclonal <- vector(mode = "numeric", length = nrow(jsi))
primary_private_subclonal <- vector(mode = "numeric", length = nrow(jsi))
shared_subclonal <- vector(mode = "numeric", length = nrow(jsi))
shared_clonal <- vector(mode = "numeric", length = nrow(jsi))

### Loop through sample pairs to calculate number of private/shared clonal/subclonal mutations
for (row in 1:nrow(jsi)) {
  # Select mutations for each sample
  tmp1 <- df %>% dplyr::filter(Tumor.ID == jsi[row, "sample1"]$sample1)
  tmp2 <- df %>% dplyr::filter(Tumor.ID == jsi[row, "sample2"]$sample2)
  
  ### Count number of private clonal mutations in metastasis
  metastasis_private_clonal_priv <- tmp2 %>% dplyr::filter(CLS %in% clonal & !(MutID %in% tmp1$MutID))
  metastasis_private_clonal[row] <- nrow(metastasis_private_clonal_priv)
  
  ### Count number of private clonal mutations in primary
  primary_private_clonal_priv <- tmp1 %>% dplyr::filter(CLS %in% clonal & !(MutID %in% tmp2$MutID))
  primary_private_clonal[row] <- nrow(primary_private_clonal_priv)
  
  ### Count number of private subclonal mutations in primary
  metastasis_private_subclonal_priv <- tmp2 %>% dplyr::filter(CLS %in% subclonal & !(MutID %in% tmp1$MutID))
  metastasis_private_subclonal[row] <- nrow(metastasis_private_subclonal_priv)
  
  ### Count number of private subclonal mutations in primary
  primary_private_subclonal_priv <- tmp1 %>% dplyr::filter(CLS %in% subclonal & !(MutID %in% tmp2$MutID))
  primary_private_subclonal[row] <- nrow(primary_private_subclonal_priv)
  
  ### Count number of shared subclonal mutation in paired samples
  shared_mutations <- intersect(tmp1$MutID, tmp2$MutID)
  shared_subclonal_priv <- tmp1 %>% dplyr::filter(CLS %in% subclonal & MutID %in% shared_mutations)
  shared_clonal_priv <- tmp1 %>% dplyr::filter(CLS %in% clonal & MutID %in% shared_mutations)
  
  shared_subclonal[row] <- nrow(shared_subclonal_priv)
  shared_clonal[row] <- nrow(shared_clonal_priv)
}

### Add data to dataframe
jsi$primary_private_clonal <- primary_private_clonal
jsi$metastasis_private_clonal <- metastasis_private_clonal
jsi$metastasis_private_subclonal <- metastasis_private_subclonal
jsi$primary_private_subclonal <- primary_private_subclonal
jsi$shared_subclonal <- shared_subclonal
jsi$shared_clonal <- shared_clonal

### Calculate Jaccard similarity index
jsi <- jsi %>% rowwise() %>%
  mutate(JSI = shared_subclonal / (metastasis_private_clonal + primary_private_clonal + shared_subclonal)) %>%
  mutate(Seeding = case_when(JSI >= 0.3 ~ "Polyclonal", TRUE ~ "Monoclonal"))

jsi.fil <- jsi %>% filter(type == "primary_metastasis")

plt <- ggplot(jsi.fil, aes(x = type, y=JSI)) + 
  geom_beeswarm(size = 3, cex = 3, color = "#2E86AB") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "black") +
  xlab("") + ylim(0, 0.35) + theme_minimal() + 
  theme(axis.text.y = element_text(size = 8), 
        axis.title = element_blank(), 
        axis.text.x = element_blank())

pdf("results/figures/mutational_timing/jaccard_similarity_index.pdf", width = 5, height = 5)
print(plt)
dev.off()

### Save results
readr::write_delim(jsi, "data/processed/mutation_timing/jaccard_similarity/mPPGL_jsi_estimates.txt", delim = "\t")
