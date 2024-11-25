#Required to run the following script before continuing: 1_data_merge_and_filter.R

# Load the required libraries
library(tidyverse)
library(ggpubr)
library(stringr)

# Calculate coverage per site
cov_plot <- combined_data_filtered %>%
  ggplot(aes(x = POS_ID, y = cov, color = PCR_prod)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0.5, outlier.size = 0.5, outlier.stroke = 0) +
  xlab("Position ID") + 
  ylab("Read deph (x)") +
  theme_classic(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
cov_plot

# Create table for median coverage per site for reporting in the text
median_cov_per_site <- combined_data_filtered %>%
  group_by(POS_ID) %>%
  summarise(median_cov = median(cov),
            min_cov = min(cov),
            max_cov = max(cov),
            lower_quantile_cov = quantile(cov)[2],
            upper_quantile_cov = quantile(cov)[4])

## Save cov_plot as .svg file
ggsave(
  "figures/cov_plot.svg",
  plot = cov_plot,
  scale = 1,
  width = 16,
  height = 10,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Clean up global environment
rm(cov_plot)
rm(median_cov_per_site)

# Calculate call rate: percentage of samples that have >100x coverage on a specific site
call_rate_table <- combined_data_filtered %>%
  mutate(called = case_when(cov >= 100 ~ 1,
                            cov < 100 ~ 0,
                            TRUE ~ 0)) %>%
  group_by(POS_ID, PCR_prod) %>%
  summarise(n_called = sum(called)) %>%
  mutate(call_rate = n_called / 365 * 100)

## Make plot of call rate
call_rate_plot <- call_rate_table %>%
  ggplot(aes(x = POS_ID, y = call_rate, color = PCR_prod, fill = PCR_prod)) +
  geom_hline(yintercept = 0, linewidth = 0.25) + 
  geom_col() + 
  xlab("Position ID") + 
  ylab("Call rate (%)") +
  theme_classic(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
call_rate_plot

## Save call_rate_plot
ggsave(
  "figures/call_rate_plot.svg",
  plot = call_rate_plot,
  scale = 1,
  width = 16,
  height = 10,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)
rm(call_rate_table)
rm(call_rate_plot)

# Load ONT data table
ONT_data_available <- read_delim("data/ONT_data_available_anonymized.txt", delim = "\t")
ONT_data_available_blood <- ONT_data_available %>%
  filter(blood_ONT_data_SMN_locus == "Yes")

# Get metrics of the dataset:
## total number of patients in ONT dataset, returns 29
length(unique(ONT_data_available$Anonymized_ID))
## total number of patients in bisulfite data, returns 365
length(unique(combined_data_filtered$Anonymized_ID))
## number of patients that overlap between ONT and bisulfite, returns 24
length(intersect(ONT_data_available$Anonymized_ID, combined_data_filtered$Anonymized_ID))
## patients in ONT dataset but not in bisulfite dataset, returns 5
length(setdiff(ONT_data_available$Anonymized_ID, combined_data_filtered$Anonymized_ID))
## total number of patients when combining ONT and bisulfite data, returns 370
length(union(ONT_data_available$Anonymized_ID, combined_data_filtered$Anonymized_ID))
## number of blood samples that overlap between ONT and bisulfite, returns 9
length(intersect(ONT_data_available_blood$Anonymized_ID, combined_data_filtered$Anonymized_ID))

# Check similarity between ONT and bisulfite data
## Reform bisulfite data
bisulfite_methylation_blood <- combined_data_filtered %>%
  select(Anonymized_ID, start, percentage) %>%
  filter(Anonymized_ID %in% intersect(ONT_data_available_blood$Anonymized_ID, combined_data$Anonymized_ID)) %>%
  mutate(technique = "bisulfite")
## Load in ONT methylation percentages for blood samples
ONT_bed_full <- read_delim("data/modbam2bed_ONT_per_patient_anonymized.bed", delim = "\t", col_names = FALSE) #load full modbam2bed output

## Change column names of ONT_bed_full
colnames(ONT_bed_full) <- c("CHROM", "start", "POS", "modbase", "score", "strand",
                        "remove1", "remove2", "remove3",
                        "coverage", "modpercentage", "Ncanonical", "Nmodified",
                        "Nfilter", "sample")

## Filter to retain blood samples only and adjust the coordinate of the '-' strand
ONT_bed_full_blood <- ONT_bed_full %>%
  select(CHROM, start, POS, strand, coverage, modpercentage, Ncanonical, Nmodified, sample) %>%
  filter(grepl("blood", sample)) %>% 
  mutate(POS = case_when(strand == "-" ~ POS - 1,
                         strand == "+" ~ POS)) %>% #change coordinate of the - strand rows, so that they correspond to the + strand and can be merged later
  mutate(Anonymized_ID = str_extract(sample, "SMA_\\d{2,3}")) %>%
  mutate(start = POS)

rm(ONT_bed_full)

## Merge the data from + and - strands to get 1 value for methylation percentage per CpG site
ONT_bed_strands_merged_blood <- ONT_bed_full_blood %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "chr"))) %>%
  group_by(Anonymized_ID, CHROM, start, POS) %>%
  summarise(Nmod = sum(Nmodified),
            Ncan = sum(Ncanonical),
            cov_tot = sum(coverage)) %>%
  mutate(cov_called = Ncan + Nmod,
         percentage = Nmod / cov_called * 100) #attention: this percentage is calculated over all bases that have a methylation call. If you want to include the filtered bases in the total, divide by cov_tot
rm(ONT_bed_full_blood)

# Ungroup dataframe and filter for Anonymized_ID in both dataframes available blood and combined_data
ONT_methylation_blood <- ONT_bed_strands_merged_blood %>%
  ungroup() %>%
  select(Anonymized_ID, start, percentage) %>%
  filter(start %in% bisulfite_methylation_blood$start) %>%
  filter(Anonymized_ID %in% intersect(ONT_data_available_blood$Anonymized_ID, combined_data$Anonymized_ID)) %>%
  mutate(technique = "ONT")

rm(ONT_bed_strands_merged_blood)
# Rbind the two tables of ONT and bisulfite --> both blood
ONT_bisulfite_comparison <- rbind(ONT_methylation_blood, bisulfite_methylation_blood)

# Compare the two methods in a correlation plot
ONT_bisulfite_comparison_spread <- ONT_bisulfite_comparison %>%
  spread(technique, percentage) %>%
  filter(!is.na(ONT) & !is.na(bisulfite))

ONT_bisulfite_correlation <- ONT_bisulfite_comparison_spread %>%
  ggplot(aes(x = ONT, y = bisulfite)) + 
  geom_point(size = 1, alpha = 0.5, stroke = 0) +
  stat_smooth(colour = "black", method = "lm", linewidth = 0.5) +
  xlab("ONT methylation percentage (%)") +
  ylab("Bisulfite methylation percentage (%)") +
  theme_classic(base_size = 7.5) +
  theme(legend.position = "none")
ONT_bisulfite_correlation

ggsave(
  "figures/ONT_bisulfite_correlation.svg",
  plot = ONT_bisulfite_correlation,
  scale = 1,
  width = 5,
  height = 5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Statistical correlation test
## Check assumptions for pearson correlation first
### 1. is the covariation linear? -> yes
### 2. does the data of each of the 2 variables follow a normal distribution?
shapiro.test(ONT_bisulfite_comparison_spread$ONT)
shapiro.test(ONT_bisulfite_comparison_spread$bisulfite)
### 3. The data is not normally distributed, so use spearman correlation test
cor.test(ONT_bisulfite_comparison_spread$ONT, ONT_bisulfite_comparison_spread$bisulfite, method="spearman")

#clean up variables
rm(ONT_bisulfite_comparison)
rm(ONT_bisulfite_comparison_spread)
rm(ONT_bisulfite_correlation)
rm(ONT_data_available)
rm(ONT_data_available_blood)
rm(ONT_methylation_blood)
rm(bisulfite_methylation_blood)
