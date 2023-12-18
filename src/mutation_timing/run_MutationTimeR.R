library(VariantAnnotation)
library(tidyverse)
library(MutationTimeR)

### Load SNV and CNV data
snv <- readr::read_csv("data/raw/snvs/PP002-DZ3A.ffpe_scenario.local-fdr.paired_somatic.vep.report.csv")

### Create fields required for rowRanges object
snv <- snv %>% dplyr::mutate(VCFStrand = case_when(STRAND == 1 ~ "+", STRAND == -1 ~ "-")) %>%
  dplyr::mutate(VCFString = stringr::str_c(Chr, Start, VCFStrand, sep = ":"))

### Create GRanges object used for rowRanges object
snv_rowRanges <- GRanges(snv$VCFString, REF = DNAStringSet(snv$REF),
                         ALT = snv$ALT, QUAL = NA, FILTER = snv$FILTER)

### Create info portion of GRanges object
snv_counts <- snv %>% dplyr::select(Tumor.Depth, Tumor.AltDepth) %>%
  dplyr::rename(t_ref_count = Tumor.Depth, t_alt_count = Tumor.AltDepth)

### Create SNV input object for MutationTimeR
mutation_time_snv_input <- VariantAnnotation::expand(VCF(rowRanges = snv_rowRanges, info = DataFrame(snv_counts)))

### Load CNV data
cnv <- readr::read_delim("data/raw/cnvs/sequenza/PP002-DZ3A_segments.txt", delim = "\t")

cnv <- cnv %>% dplyr::mutate(ranges = stringr::str_c(start.pos, end.pos, sep = "-")) %>%
  dplyr::mutate(VCFString = stringr::str_c(chromosome, ranges, "*", sep = ":"))

mutation_time_cnv_input <- GRanges(cnv$VCFString, major_cn = cnv$A, 
                    minor_cn = cnv$B, clonal_frequency = cnv$Bf)

mt <- mutationTime(mutation_time_snv_input, mutation_time_cnv_input, n.boot=200)

V <- as.data.frame(mt$V) %>% dplyr::select(!pAllSubclones) %>%
  dplyr::mutate(CNID = as.numeric(CNID)) %>%
  dplyr::mutate(CN = MajCN + MinCN)

T <- as.data.frame(mt$T[-c(1)])

# Add mutation timing info for plotting
mutation_time_snv_input <- addMutTime(mutation_time_snv_input, mt$V)
mcols(mutation_time_cnv_input) <- cbind(mcols(mutation_time_cnv_input),mt$T)

# Plot
pdf(file="results/mutation_timing/MutationTimeR/PP002-DZ3A-MutationTimeR.pdf", height=8.5, width=11)
plotSample(mutation_time_snv_input,mutation_time_cnv_input)
dev.off()

# Save timing and variant information
write.csv(T, 'data/processed/mutation_timing/MutationTimeR/PP002-DZ3A_mT.csv') # Save timing information
write.csv(V, 'data/processed/mutation_timing/MutationTimeR/PP002-DZ3A_mV.csv') # Save variant information

# Join all data frames
snv <- cbind(snv, V)

cnv$CNID <- 1:nrow(cnv)
snv <- dplyr::left_join(snv, cnv, by="CNID")

write.csv(snv, "data/processed/mutation_timing/MutationTimeR/PP002-DZ3A_MutTimeR.csv")