## ---------------------------
## Script Name: calculate_hrd_stats.R
## Description: Calculate HRD score for mPPHL
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-26
## ---------------------------

# Load packages
library(tidyverse)
library(HRDex)
library(ggpubr)

# Define variables
sdh <- c("SDHA", "SDHB", "SDHC")

# Load data ---------------------------------------------------------------

cnv.files <- list.files(
  "data/raw/cnvs/annotSV/", pattern = "*.csv", full.names = FALSE, recursive = FALSE
  )

# Read all CSV files into a list
df <- map(cnv.files, function(file) {
  result <- read_csv(paste0("data/raw/cnvs/annotSV/", file))
  return(result)
}) %>% bind_rows()

calculate_HRD_score <- function(df) {
  
  # Extract sample name from data
  sample <- unique(df$TumorID)
  
  df <- df %>% 
    select(SV.Chrom, SV.Start, SV.End, Segment.CN, Segment.CNa, Segment.CNb) %>%
    distinct(SV.Chrom, SV.Start, SV.End, Segment.CN, Segment.CNa, Segment.CNb) %>%
    rename(chromosome = SV.Chrom, start.pos = SV.Start, end.pos = SV.End,
           CNt = Segment.CN, A = Segment.CNa, B = Segment.CNb) %>%
    mutate(chromosome = paste0("chr", chromosome))
  
  # Get seq and CN dataframes
  seq.dat <- HRDex::preprocessHRD(df, "grch38")
  CN.dat <- HRDex::getCNt(seq.dat)
  
  # Compute the scores
  HRD.Score.sum <- getHRD.Score(seq.dat, CN.dat, type = "sum")
  HRD.Score.average <- getHRD.Score(seq.dat, CN.dat, type = "average")
  LST.raw <- getLST(seq.dat)
  LOH.raw <- getLOH(seq.dat)
  NTAI.norm <- getNTAI.norm(seq.dat, CN.dat)
  NTAI.raw <- getNTAI.raw(seq.dat)
  
  # Put results together
  hrd.stats <- data.frame(
    Sample = sample,
    HRD.Score.sum = HRD.Score.sum,
    HRD.Score.average = HRD.Score.average,
    LST.raw = LST.raw,
    LOH.raw = LOH.raw,
    NTAI.norm = NTAI.norm,
    NTAI.raw = NTAI.raw
  )
  
  return(hrd.stats)

}

hrd.stats <- df %>%
  group_by(TumorID) %>%
  group_split() %>%
  map_df(calculate_HRD_score)

# SAVE RESULTS
write_csv(hrd.stats, "data/processed/cnvs/HRDex/mPPGL.hrd.scores.csv")

# Plot HRD Statistics -----------------------------------------------------

# Load metadata and merge
meta <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")
hrd.stats <- left_join(hrd.stats, meta, by=c("Sample" = "sample"))

# Create first plot comparing primary and metastatic tumors
plt1<- ggboxplot(data = hrd.stats, x="category", y="HRD.Score.average", fill="category",
                          palette = c('#083D77', "#F4D35E"), 
                          ylab = "Mean HRD Score", xlab="", add="jitter", add.params = list(color="#071E22", size=2)) +
  scale_x_discrete(labels = c("Primary", "Metastatic"))


plt1 <- plt1 + stat_compare_means(method = "t.test", comparisons = list(c("Primary", "Metastatic")))

plt1 <- plt1 + theme(legend.position = "none", 
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.border = element_blank(),
                     axis.text.x = element_text(size = 8, color = "black"),
                     axis.text.y = element_text(size = 8, color = "black"),
                     axis.title.x = element_text(size = 8),
                     axis.title.y = element_text(size = 8),
                     axis.line.y = element_line())


# Create plot comparing HRD scores by germline mutations
hrd.stats <- hrd.stats %>% 
  filter(germline != "RET")

### MORE STATS
aov1 <- aov(HRD.Score.average ~ germline, data=hrd.stats)
summary(aov1)

plt2 <- ggboxplot(data = hrd.stats, x="germline", y="HRD.Score.average", fill="germline",
                           order = c("WT", "SDHA", "SDHB", "SDHC"), palette = c("#D11149", "#1A8FE3", "#F17105", "#CECCCC"),
                           ylab = "", xlab="", add="jitter", add.params = list(color="#071E22", size=2))


plt2 <- plt2 + stat_compare_means(method = "wilcox.test", comparisons = list(c("WT", "SDHA"), 
                                                                               c("WT", "SDHB"),
                                                                               c("WT", "SDHC"),
                                                                               c("SDHA", "SDHB"),
                                                                               c("SDHA", "SDHC"),
                                                                               c("SDHB", "SDHC")))

plt2 <- plt2 + theme(legend.position = "none", 
                       strip.background = element_blank(),
                       strip.text = element_blank(),
                       panel.border = element_blank(),
                       axis.text.x = element_text(size = 8, color = "black"),
                       axis.text.y = element_text(size = 8, color = "black"),
                       axis.title.x = element_text(size = 8),
                       axis.title.y = element_text(size = 8),
                       axis.line.y = element_line())

fig1 <- cowplot::plot_grid(plt1, plt2, rel_widths = c(0.8, 1), labels = "auto", label_size = 12)

pdf("results/hrd_scores/hrd_scores_summary.pdf", width=7.5, height=5.625)
print(fig1)
dev.off()
