## ---------------------------
## Script Name: calculate_arm_level_changes.R
## Description: Calculate the percentage of an arm affected by gain, amplification,
## deletion, loss, or neutral copy number event
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-09-26
## ---------------------------

library(tidyverse)
library(optparse)

# PARSE COMMAND LINE OPTIONS
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "Path to data frame",
              metavar = "character"
  ),
  make_option(c("-c", "--cytobands"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-p", "--ploidy"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-o", "--output"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  ),
  make_option(c("-s", "--summary"),
              type = "character",
              default = "out.txt",
              help = "output file name [default = %default]",
              metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

df <- readr::read_delim(opt$input, delim = "\t")
cytos <- readr::read_delim(opt$cytobands)
ploidy <- readxl::read_xlsx(opt$ploidy, sheet = "tumors")

### SELECT UNIQUE EVENTS
df <- df %>% dplyr:::mutate(Segment.Type = factor(Segment.Type, levels = c("del", "loss", "amp", "gain", "neutral"))) %>%
  dplyr::distinct(TumorID, AnnotSV.ID, SV.Chrom, SV.Start, SV.End, Segment.CN, Segment.Length, Segment.Zyg, Segment.Type) %>%
  dplyr::mutate(Segment.Length = as.numeric(stringr::str_remove(Segment.Length, "bp"))) %>%
  dplyr::full_join(cytos, by = c("SV.Chrom" = "Chrom"), relationship = "many-to-many") 

### MERGE CYTOBAND INFORMATION AND SELECT APPROPRIATE ARM
df <- df %>% dplyr::filter(SV.Start >= Start & SV.Start < End) %>%
  dplyr::mutate(Arm = factor(Arm, levels = cytos$Arm)) %>% 
  dplyr::left_join(ploidy, by = c("TumorID"="sample"))

### CALCULATE PERCENT OF ARM LOSS
df <- df %>%
  dplyr::mutate(Percent.Change = case_when(
    Segment.Type %in% c("amp", "gain") & Segment.CN <= ploidyRounded ~ 0,
    Segment.Type == "neutral" & Segment.Zyg != "loh" ~ 0,
    TRUE ~ (Segment.Length / Arm.Length) * 100
  )) %>%
  dplyr::group_by(TumorID, Arm, Segment.Type, .drop = FALSE) %>%
  dplyr::summarise(Total = sum(Percent.Change)) %>%
  reshape2::dcast(formula = TumorID + Arm ~ Segment.Type, fill = 0)

### CALCULATE BINARY GAIN OR LOSS VARIABLE
arm_changes <- df %>% dplyr::mutate(Arm.Gain = case_when(amp >= 50 ~ 1, gain >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(Arm.Loss = case_when(del >= 50 ~ 1, loss >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(LOH = case_when( neutral >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::select(TumorID, Arm, Arm.Gain, Arm.Loss) %>%
  dplyr::rename(Sample = TumorID)

### SAVE RESULTS
readr::write_delim(df, opt$output, delim="\t")
readr::write_delim(arm_changes, opt$summary, delim="\t")