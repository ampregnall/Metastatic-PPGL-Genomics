library(tidyverse)
library(viridis)

### Load data
df <- readr::read_delim("data/processed/mutation_timing/ccf_estimates/mPPGL_cancer_cell_fractions.txt")

### Define plotting cohorts
group1 <- sort(unique(df$Tumor.ID))[1:12]
group2 <- sort(unique(df$Tumor.ID))[13:24]
group3 <- sort(unique(df$Tumor.ID))[25:36]
group4 <- sort(unique(df$Tumor.ID))[37:48]

### Plot distribution of CCF estimates in overall cohort
plot <- ggplot(df, aes(x=CCF)) + geom_density(fill="lightblue") +
  xlab("Cancer Cell Fraction") + ylab("Density") + labs(title="Distribution of CCF Estimates in mPPGL Cohort") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "#D36492") +
  theme_minimal() + theme(legend.position = "none", 
       plot.title = element_text(size = 16, hjust=0.5),
       axis.text.x = element_text(size = 16),
       axis.text.y = element_text(size = 16),
       axis.title.x = element_text(size = 16),
       axis.title.y = element_text(size = 16))

### Save results
pdf("results/figures/ccf_analysis/ccf_distribution.pdf", width = 16, height = 8.5)
print(plot)
dev.off()


# Distribution of CCF Estimates by Sample ---------------------------------

### Subset data for plotting on multiple pages
df1 <- df %>% dplyr::filter(Tumor.ID %in% group1)
df2 <- df %>% dplyr::filter(Tumor.ID %in% group2)
df3 <- df %>% dplyr::filter(Tumor.ID %in% group3)
df4 <- df %>% dplyr::filter(Tumor.ID %in% group4)

### Create plot for first page
plot1 <- ggplot(df1, aes(x=CCF)) + geom_density(fill="lightblue") +
  xlab("Cancer Cell Fraction") + ylab("Density") + facet_wrap(~Tumor.ID) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "#D36492") +
  theme_minimal() + theme(legend.position = "none", 
                          plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 16),
                          axis.title.y = element_text(size = 16),
                          strip.text = element_text(size = 16))

### Create plot for second page
plot2 <- ggplot(df2, aes(x=CCF)) + geom_density(fill="lightblue") +
  xlab("Cancer Cell Fraction") + ylab("Density") + facet_wrap(~Tumor.ID) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "#D36492") +
  theme_minimal() + theme(legend.position = "none", 
                          plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 16),
                          axis.title.y = element_text(size = 16),
                          strip.text = element_text(size = 16))

### Create plot for third page
plot3 <- ggplot(df3, aes(x=CCF)) + geom_density(fill="lightblue") +
  xlab("Cancer Cell Fraction") + ylab("Density") + facet_wrap(~Tumor.ID) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "#D36492") +
  theme_minimal() + theme(legend.position = "none", 
                          plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 16),
                          axis.title.y = element_text(size = 16),
                          strip.text = element_text(size = 16))

### Create plot for fourth page
plot4 <- ggplot(df4, aes(x=CCF)) + geom_density(fill="lightblue") +
  xlab("Cancer Cell Fraction") + ylab("Density") + facet_wrap(~Tumor.ID) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "#D36492") +
  theme_minimal() + theme(legend.position = "none", 
                          plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 16),
                          axis.title.y = element_text(size = 16),
                          strip.text = element_text(size = 16))

### Save results
pdf("results/figures/ccf_analysis/ccf_distribution_by_sample.pdf", width = 16, height = 8.5, onefile = TRUE)
print(plot1)
print(plot2)
print(plot3)
print(plot4)
dev.off()

# Plot Proportion of Clonal/Subclonal Mutation by Sample ------------------

### Classify mutations as clonal/subclonal based on CCF
df <- df %>% dplyr::mutate(CLS = case_when(CCF.95 == 1 ~ "Clonal", TRUE ~ "Subclonal"))

### Count number of clonal/subclonal mutations in cohort
CLS_summary <- df %>% dplyr::group_by(CLS) %>%
  dplyr::summarise(count = n())

### Count number of clonal/subclonal mutations by sample
CLS_sample_summary <- df %>% dplyr::group_by(Tumor.ID, CLS) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = CLS, values_from = count) %>%
  tidyr::replace_na(list(Subclonal = 0))

# Convert data to long format for ggplot2
df_long <- CLS_sample_summary %>% gather(key = "type", value = "count", -Tumor.ID)

# Calculate proportions
df_long <- df_long %>% group_by(Tumor.ID) %>%
  mutate(proportion = count / sum(count))

# Create the bar plot
plot <- ggplot(df_long, aes(x = Tumor.ID, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of Clonal and Subclonal Mutations by Sample",
       x = "Sample", y = "Proportion", fill = "Type") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("lightblue", "#D36492")) +
  theme_minimal() + theme(plot.title = element_text(size = 16, hjust=0.5),
                          axis.text.x = element_text(size = 12, angle = 90),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 16),
                          axis.title.y = element_text(size = 16),
                          legend.title = element_text(size = 16),
                          legend.text = element_text(size = 12))

### Save results
pdf("results/figures/ccf_analysis/proportion_clonal_mutations.pdf", width = 16, height = 8.5)
print(plot)
dev.off()

# Relationship Between Tumor.AltDepth and CCF -----------------------------

### Create categorical encoding of tumor copy number
df <- df %>% dplyr::mutate(CNt_cat = case_when(Segment.CN == 0 ~ "0", Segment.CN == 1 ~ "1", Segment.CN == 2 ~"2", 
                                               Segment.CN == 3 ~ "3", Segment.CN == 4 ~ "4", TRUE ~ "5+"))

### Plot relationship between CCF and alternative allele fraction, color coded by tumor copy number 
plot <- ggplot(df, aes(x=Tumor.AltFrac, y=CCF, color=CNt_cat)) + geom_point() + theme_minimal() +
  xlab("Alternative Allele Fraction") + ylab("Cancer Cell Fraction") +
  labs(title="Correlation Between Alternative Allele Fraction and Cancer Cell Fraction") +
  scale_color_viridis(name = "Copy Number", discrete=TRUE) +
  theme(plot.title = element_text(size = 16, hjust=0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

### Save results
pdf("results/figures/ccf_analysis/ccf_alt_fraction_scatterplot.pdf", width = 16, height = 8.5)
print(plot)
dev.off()

### Subset data for plotting on multiple pages
df1 <- df %>% dplyr::filter(Tumor.ID %in% group1)
df2 <- df %>% dplyr::filter(Tumor.ID %in% group2)
df3 <- df %>% dplyr::filter(Tumor.ID %in% group3)
df4 <- df %>% dplyr::filter(Tumor.ID %in% group4)

### Create first plot
plot1 <- ggplot(df1, aes(x=Tumor.AltFrac, y=CCF, color=CNt_cat)) + geom_point() + theme_minimal() +
  xlab("Alternative Allele Fraction") + ylab("Cancer Cell Fraction") + facet_wrap(~Tumor.ID) +
  labs(title="Correlation Between Alternative Allele Fraction and Cancer Cell Fraction") +
  scale_color_viridis(name = "Copy Number", discrete=TRUE) +
  theme(plot.title = element_text(size = 16, hjust=0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

### Create second plot
plot2 <- ggplot(df2, aes(x=Tumor.AltFrac, y=CCF, color=CNt_cat)) + geom_point() + theme_minimal() +
  xlab("Alternative Allele Fraction") + ylab("Cancer Cell Fraction") + facet_wrap(~Tumor.ID) +
  labs(title="Correlation Between Alternative Allele Fraction and Cancer Cell Fraction") +
  scale_color_viridis(name = "Copy Number", discrete=TRUE) +
  theme(plot.title = element_text(size = 16, hjust=0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

plot3 <- ggplot(df3, aes(x=Tumor.AltFrac, y=CCF, color=CNt_cat)) + geom_point() + theme_minimal() +
  xlab("Alternative Allele Fraction") + ylab("Cancer Cell Fraction") + facet_wrap(~Tumor.ID) +
  labs(title="Correlation Between Alternative Allele Fraction and Cancer Cell Fraction") +
  scale_color_viridis(name = "Copy Number", discrete=TRUE) +
  theme(plot.title = element_text(size = 16, hjust=0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

plot4 <- ggplot(df4, aes(x=Tumor.AltFrac, y=CCF, color=CNt_cat)) + geom_point() + theme_minimal() +
  xlab("Alternative Allele Fraction") + ylab("Cancer Cell Fraction") + facet_wrap(~Tumor.ID) +
  labs(title="Correlation Between Alternative Allele Fraction and Cancer Cell Fraction") +
  scale_color_viridis(name = "Copy Number", discrete=TRUE) +
  theme(plot.title = element_text(size = 16, hjust=0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

### Save results
pdf("results/figures/ccf_analysis/ccf_alt_fraction_scatterplot_by_sample.pdf", width = 16, height = 8.5, onefile = TRUE)
print(plot1)
print(plot2)
print(plot3)
print(plot4)
dev.off()

