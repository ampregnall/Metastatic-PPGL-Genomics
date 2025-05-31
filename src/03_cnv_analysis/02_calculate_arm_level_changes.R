## ---------------------------
## Script Name: 02_calculate_arm_level_changes.R
## Description: Calculate the percentage of an arm affected by gain, del, or LOH CNV event
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-09-26
## ---------------------------

library(tidyverse)
library(optparse)

df <- readr::read_delim(opt$input, delim = "\t")
cytos <- readr::read_delim(opt$cytobands)
ploidy <- readxl::read_xlsx(opt$ploidy, sheet = "tumors")

arms <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q",
          "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q",
          "13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q",
          "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q", "Xp", "Xq")

### SELECT UNIQUE EVENTS
df <- df %>% dplyr:::mutate(Segment.Type = factor(Segment.Type, levels = c("del", "loss", "amp", "gain", "neutral"))) %>%
  dplyr::distinct(TumorID, AnnotSV.ID, SV.Chrom, SV.Start, SV.End, Segment.CN, Segment.Length, Segment.Zyg, Segment.Type) %>%
  dplyr::mutate(Segment.Length = as.numeric(stringr::str_remove(Segment.Length, "bp"))) %>%
  dplyr::left_join(cytos, by = c("SV.Chrom" = "chrom")) 

# Handle Sequenza annotations that are split across chromosomal arms
df <- df %>%
  # Create rows for both conditions at once with case_when
  dplyr::filter(SV.Start >= p_start & SV.Start <= p_end & SV.End > q_start & SV.End <= q_end) %>% # Select events that span arms and split across arm
  dplyr::mutate(
    SV.Start = case_when(
      SV.End > q_start & SV.End <= q_end ~ q_start,
      TRUE ~ SV.Start
    ),
    SV.End = case_when(
      SV.Start >= p_start & SV.Start <= p_end ~ p_end,
      TRUE ~ SV.End
    ),
    Segment.Length = SV.End - SV.Start
  ) %>%
  # Combine with the filtered original dataset
  dplyr::bind_rows(
    df %>%
      filter((SV.Start >= p_start & SV.End <= p_end) | (SV.Start >= q_start & SV.End <= q_end))
  )

# Annotate with chromosomal arm
df <- df %>% dplyr::mutate(arm = case_when(SV.Start >= p_start & SV.End <= p_end ~ str_c(SV.Chrom, "p"),
                                SV.Start >= q_start & SV.End <= q_end ~ str_c(SV.Chrom, "q")),
                       arm_length = case_when(SV.Start >= p_start & SV.End <= p_end ~ p_length,
                                                    SV.Start >= q_start & SV.End <= q_end ~ q_length)) %>% 
  dplyr::mutate(arm = factor(arm, levels = arms)) %>% 
  dplyr::left_join(ploidy, by = c("TumorID"="sample"))

# Calculate proportion of arm affected by amp/del/loh event
df <- df %>% dplyr::mutate(Percent.Change = case_when(
    Segment.Type %in% c("amp", "gain") & Segment.CN <= ploidyRounded ~ 0,
    Segment.Type == "neutral" & Segment.Zyg != "loh" ~ 0,
    TRUE ~ (Segment.Length / arm_length) * 100
  )) %>%
  dplyr::group_by(TumorID, arm, Segment.Type, .drop = FALSE) %>%
  dplyr::summarise(Total = sum(Percent.Change)) %>%
  reshape2::dcast(formula = TumorID + arm ~ Segment.Type, fill = 0)

### CALCULATE BINARY GAIN OR LOSS VARIABLE
arm_changes <- df %>% dplyr::mutate(Arm.Gain = case_when(amp >= 50 ~ 1, gain >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(Arm.Loss = case_when(del >= 50 ~ 1, loss >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::mutate(LOH = case_when( neutral >= 50 ~ 1, TRUE ~ 0)) %>%
  dplyr::select(TumorID, arm, Arm.Gain, Arm.Loss, LOH) %>%
  dplyr::rename(Sample = TumorID)

### SAVE RESULTS
readr::write_delim(df, opt$output, delim="\t")
readr::write_delim(arm_changes, opt$summary, delim="\t")