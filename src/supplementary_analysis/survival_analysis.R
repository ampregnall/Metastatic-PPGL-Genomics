## ---------------------------
## Script Name: survival_analysis.R
## Description: Perform survival analysis
##
## Author: Andrew M. Pregnall
## Email: andrew.pregnall@pennmedicine.upenn.edu
## 
## Date Created: 2025-01-13
## ---------------------------

library(tidyverse, quietly = TRUE)
library(survival)
library(survminer)
library(lubridate)
library(readxl)

# Load data
df <- read_xlsx("metadata/final_sample_list_survival.xlsx", sheet = "Survival")

# Survival Analysis --------------------------------------------

# Overall survival
sfit <- survfit(Surv(time = os_time, event = deceased_status_coded) ~ 1, data=df)
s1 <- ggsurvplot(sfit, risk.table = TRUE, palette = "#2E86AB", tables.theme = theme_cleantable()) 

s1$plot <- s1$plot + labs(title="Overall survival") + xlab("Time (years)") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 6))

s1$table <- s1$table + xlab("Time (years)") + ylab("") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 8))

# Recurrence free survival
sfit2 <- survfit(Surv(time = metastasis_time, event = metastasis_event) ~ 1, data=df)
s2 <- ggsurvplot(sfit2, risk.table = TRUE, palette = "#2E86AB", tables.theme = theme_cleantable()) 

s2$plot <- s2$plot + labs(title="Recurrence free survival") + xlab("Time (years)") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 6))

s2$table <- s2$table + xlab("Time (years)") + ylab("") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size = 6))

fig <- cowplot::plot_grid(s2$plot, s1$plot, s2$table, s1$table, nrow = 2, 
                          rel_heights = c(0.8, 0.2), labels = c("A", "B", "", ""), align = "hv")

pdf("results/figures/survival.pdf", width = 15, height = 10)
print(fig)
dev.off()
