## ---------------------------
## Script Name: NAME
## Description: PURPOSE
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-26
## ---------------------------

# Load packages
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)


df <- readxl::read_xlsx("results/ibd_scores/PPGL-IBD.xlsx", sheet = "paired")
df$cat <- "Tumor Pairs" # Dummy variable for plotting
tumor_metadata <- readxl::read_xlsx("metadata/mPPGL-Metadata.xlsx", sheet = "tumors")

df2 <- readr::read_delim("data/processed/snvs/drivers/PPGL.somatic.data.combined.filtered.txt")
df2 <- df2 %>% dplyr::select(c(1:6, 8))
df2 <- df2 %>% dplyr::left_join(tumor_metadata, by=c("Tumor.ID"="sample")) %>% 
  dplyr::filter(cohort == "Derivation") %>% 
  dplyr::mutate(VarID = stringr::str_c(Normal.ID, Chr, Start, REF, ALT, Gene))

df3 <- df %>% 
  dplyr::filter(IID1 %in% df2$Tumor.ID) %>% 
  dplyr::select(IID1, IID2)

shared_vars <- vector(mode = "numeric", length = length(df3))
total_vars <- vector(mode = "numeric", length = length(df3))
tumor_1_vars <- vector(mode = "numeric", length = length(df3))
tumor_2_vars <- vector(mode = "numeric", length = length(df3))

for (row in 1:nrow(df3)) {
  tmp1 <- df2 %>% dplyr::filter(Tumor.ID == df3[row, "IID1"]$IID1)
  tmp2 <- df2 %>% dplyr::filter(Tumor.ID == df3[row, "IID2"]$IID2)
  tumor_1_vars[row] <- nrow(tmp1)
  tumor_2_vars[row] <- nrow(tmp2)
  shared_vars[row] <- length(intersect(tmp1$VarID, tmp2$VarID))
  total_vars[row] <- length(tmp1$VarID) + length(tmp2$VarID) - shared_vars
}

df3$shared_vars <- shared_vars
df3$total_vars <- total_vars
df3$tumor_1_vars <- tumor_1_vars
df3$percent <- df3$shared_vars / df3$tumor_1_vars * 100
df3$cat <- "Tumor Pairs"

median(df3$percent)
summary(df3$percent)

plt <- ggplot(df3, aes(x = cat, y=percent)) + 
  geom_beeswarm(size = 3, cex = 3, color = "#7260a6") + 
  geom_beeswarm(data = df, aes(x = cat, y = PI_HAT * 100), cex = 2, size = 3, color = "#6a9e3b") +
  labs(title="Relatedness", y = "") +
  theme_minimal() + theme(strip.background = element_blank(),
                          strip.text = element_blank(),
                          panel.border = element_blank(),
                               plot.title = element_text(size = 8, hjust = 0.5),
                               axis.title.x = element_text(size = 6, color = "black"),
                               axis.text.y = element_text(size = 6, color = "black"),
                               axis.title.y = element_blank(),
                               legend.text = element_text(size = 16),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(), 
                               legend.position = "none",
                               plot.margin = unit(c(0, 0, 0, 0), "cm"))

pdf("results/figures/clonality_analysis/relatedness.pdf", width = 1.5, height = 1.5)
print(plt)
dev.off()