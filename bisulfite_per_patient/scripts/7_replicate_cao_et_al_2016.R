# Run the following script before continuing: 1_data_merge_and_filter.R
# This script checks the methylation status of the sites found in Cao et al, 2016 (PMC4710843)

# Load required libraries
library(rstatix)
library(ggpubr)

# Define sites from Cao et al, 2016 (PMC4710843)
cao <- combined_data_filtered %>%
  mutate(Cao_site = case_when(stop == 71381014 ~ "CGI1_02",
                              stop == 71381150 ~ "CGI1_08",
                              stop == 71381602 ~ "CGI2_09", #CGI2_09 includes 3 CpG sites
                              stop == 71381604 ~ "CGI2_09", #CGI2_09 includes 3 CpG sites
                              stop == 71381607 ~ "CGI2_09", #CGI2_09 includes 3 CpG sites
                              stop == 71382781 ~ "CGI4_02",
                              stop == 71382829 ~ "CGI4_03",
                              stop == 71382879 ~ "CGI4_06",
                              stop == 71382890 ~ "CGI4_07",
                              stop == 71382955 ~ "CGI4_10"
                              )) %>%
  filter(!is.na(Cao_site)) %>%
  filter(SMN2_CN == 3) %>%
  filter(SMN1_CN == 0) %>%
  filter(mut_859_GC == "no") %>%
  filter(SMA_type_numeric == 1 | SMA_type_numeric == 2 | SMA_type_numeric == 3)

# number of patients with type 1
length(unique(filter(cao, SMA_type_numeric == 1)$Anonymized_ID))
# number of patients with type 2
length(unique(filter(cao, SMA_type_numeric == 2)$Anonymized_ID))
# number of patients with type 3
length(unique(filter(cao, SMA_type_numeric == 3)$Anonymized_ID))

cao_plot <- cao %>% ggplot(aes(x = Cao_site, y = percentage, color = as.factor(SMA_type_numeric))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), alpha = 0.2, size = 0.2) +
  xlab("CpG site") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(name = "SMA type", labels = c("Type 1", "Type 2", "Type 3"),
                     values  = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_classic(base_size = 7.5)

cao_plot

# Save cao_plot as .svg
ggsave(
  "figures/cao_plot.svg",
  plot = cao_plot,
  scale = 1,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Statistics
## Check assumptions for one-way ANOVA for 1 site
cao_CGI1_02 <- cao %>%
  filter(Cao_site == "CGI1_02")
## Identify outliers
cao_CGI1_02 %>%
  group_by(SMA_type_numeric) %>%
  identify_outliers(percentage)
## determine if it is normally distributed
cao_CGI1_02 %>%
  group_by(SMA_type_numeric) %>%
  shapiro_test(percentage)
ggqqplot(cao_CGI1_02, x = "percentage", facet.by = "SMA_type_numeric")

## data is not normally distributed --> do kruskal-wallis test for all sites, do multiple testing correction because we are testing 8 sites
kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI1_02"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI1_02"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI1_08"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI1_08"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI2_09"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI2_09"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_02"))# p = 0.035 without multiple testing correction
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_02"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_03"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_03"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_06"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_06"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_07"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_07"))$p.value, method = "fdr", n = 8)

kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_10"))
p.adjust(kruskal.test(percentage ~ as.factor(SMA_type_numeric), data = filter(cao, Cao_site == "CGI4_10"))$p.value, method = "fdr", n = 8)

