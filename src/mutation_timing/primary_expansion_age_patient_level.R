## ---------------------------
## Script Name: primary_tumor_expansion_age.R
## Description: Estimate how long before diagnosis pPCC and pPGL originated
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2024-10-12
## ---------------------------

library(tidyverse)
library(spatstat)

### Load data on doubling time and diameter
df <- readxl::read_xlsx("metadata/PPGL_doubling_time_estimates.xlsx", sheet = "tumors")

limit_size <- 1e11
log_limit_size <- log(limit_size)

df <- df %>%
  mutate(
    dx_size = tumor_size_t1 * 1e8,
    doubling_time_days = doubling_time_yr * 365,
    # Calculate beta (growth rate)
    beta = 1/doubling_time_days * log((log_limit_size - log(dx_size)) / (log_limit_size - log(2 * dx_size))),
    # Calculate alpha (growth parameter)
    alpha = log_limit_size * beta,
    # Calculate expected age in years
    exp_age = -1/beta * log(1 - log(dx_size) * beta / alpha) / 365
  )

# Filter to non-negative estimates and those less than patient age
df.fil <- df %>% filter(exp_age > 0 & exp_age < age)
