library(VariantAnnotation)
library(tidyverse)
library(MutationTimeR)

# Load all SNV files in directory
snv.files <- list.files(
  path = "data/raw/snvs/", pattern = "*.csv", full.names = FALSE, recursive = FALSE
)

# Load all CNV files in directory
cnv.files <- list.files(
  path = "data/raw/cnvs/sequenza/", pattern = "*.txt", full.names = FALSE, recursive = FALSE
)

# Load sample metadata
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

# Get samples names from both lists
snv.samples <- sapply(snv.files, str_extract, "[^.]+")
cnv.samples <- sapply(cnv.files, str_extract, "[^_]+")

# Get matching samples in both lists
samples <- intersect(snv.samples, cnv.samples)
missing <- setdiff(snv.samples, cnv.samples) # Sequenza files are missing

# Remove missing samples from SNV list
snv.files <- snv.files[snv.samples %in% missing == FALSE]
snv.samples <- snv.samples[snv.samples %in% missing == FALSE]

# Ensure lists are ordered the same
snv.index <- match(samples, snv.samples)
cnv.index <- match(samples, cnv.samples)
snv.files <- snv.files[snv.index]
cnv.files <- cnv.files[cnv.index]

# Load data
snvs <- map(snv.files, function(file) {
  read_delim(paste0("data/raw/snvs/", file),
    delim = ",",
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
})

cnvs <- map(cnv.files, function(file) {
  read_delim(paste0("data/raw/cnvs/sequenza/", file),
    delim = "\t",
    col_types = "cnnnnnnnnnnnn"
  )
})

# Create SNV input --------------------------------------------------------

process.MutTimeR.snvs <- function(snv) {
  
  # Create fields required for rowRanges object
  snv <- snv %>%
    mutate(VCFStrand = case_when(STRAND == 1 ~ "+", STRAND == -1 ~ "-")) %>%
    mutate(VCFString = stringr::str_c(Chr, Start, VCFStrand, sep = ":"))

  # Create GRanges object used for rowRanges object
  rowRanges <- GRanges(snv$VCFString,
    REF = DNAStringSet(snv$REF),
    ALT = snv$ALT, QUAL = NA, FILTER = snv$FILTER
  )

  # Create info portion of GRanges object
  counts <- snv %>%
    dplyr::select(Tumor.Depth, Tumor.AltDepth) %>%
    dplyr::rename(t_ref_count = Tumor.Depth, t_alt_count = Tumor.AltDepth)

  # Create SNV input object for MutationTimeR
  return(VariantAnnotation::expand(VCF(rowRanges = rowRanges, info = DataFrame(counts))))
}

snvs.mtr <- map(snvs, process.MutTimeR.snvs)

# Create CNV input --------------------------------------------------------

process.MutTimeR.cnvs <- function(cnv, i) {
  sample <- samples[i]
  purity <- meta[meta$sample == sample,]$cellularity

  cnv <- cnv %>%
    mutate(ranges = stringr::str_c(start.pos, end.pos, sep = "-")) %>%
    mutate(VCFString = stringr::str_c(chromosome, ranges, "*", sep = ":"))

  return(GRanges(cnv$VCFString, major_cn = cnv$A, minor_cn = cnv$B, clonal_frequency = purity))
}

cnvs.mtr <- map2(cnvs, seq_along(cnvs), process.MutTimeR.cnvs)

# Run MutationTimeR -------------------------------------------------------

for (i in seq_along(samples)) {
  print(paste("Running Sample:" , samples[i]))
  
  # Extract values
  s <- snvs.mtr[[i]]
  c <- cnvs.mtr[[i]]
  gender <- meta[meta$sample == samples[i],]$gender
  duplicated <- meta[meta$sample == samples[i],]$wholeGenomeDuplication
  
  # Run MutationTimeR
  mt <- mutationTime(s, c, gender = gender, isWgd = duplicated)
  
  # Extract variant information
  V <- as.data.frame(mt$V) %>% 
    dplyr::select(!pAllSubclones) %>%
    dplyr::mutate(CNID = as.numeric(CNID)) %>%
    dplyr::mutate(CN = MajCN + MinCN)
  
  # Extract CNV timing information
  T <- as.data.frame(mt$T[-c(1)])
  
  # Save timing and variant information
  print("Saving data...")
  write.csv(T, paste0("data/processed/mutation_timing/MutationTimeR/timing/", samples[i], "_mT.csv")) 
  write.csv(V, paste0("data/processed/mutation_timing/MutationTimeR/variants/", samples[i], "_mV.csv")) 
  
  # Add mutation timing info for plotting
  s <- addMutTime(s, mt$V)
  mcols(c) <- cbind(mcols(c), mt$T)
  
  # Plot
  print("Saving plot...")
  pdf(file=paste0("data/processed/mutation_timing/MutationTimeR/plots/", samples[i], "_mutationTime_plot.pdf"), height=8.5, width=11)
  plotSample(s,c)
  dev.off()
  
  # Join all data frames
  df <- cbind(snvs[i], V)
  cnvs[[i]]$CNID <- 1:nrow(cnvs[[i]])
  df <- left_join(df, cnvs[[i]], by="CNID")
  
  # Save data
  write.csv(df, paste0("data/processed/mutation_timing/MutationTimeR/variant_timing/", samples[i], "_MutTimeR.csv"))
  
  }
  