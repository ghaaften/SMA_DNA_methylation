# Run the following script before continuing: 1_data_merge_and_filter.R

# Load required libraries for "2_ONT_baseline_characteristics.R"
library(tidyverse)
library(table1)

# Create baseline table from the bed_strands_merged dataframe (script: 1_data_merge_and_filter.R)
baseline_table <- bed_strands_merged %>%
  select(SMA_ID, SMA_ID_tissue, tissue, Sex, SMN2_CN, SMA_type, SMA_type_ext, age_at_sampling, age_at_onset_years) %>%
  distinct() %>%
  mutate(Total = 1)

# Adjust columns to factors with correct labels
## Column: Sex
baseline_table$Sex <- 
  factor(baseline_table$Sex, levels=c("M","F"),
         labels=c("Male", 
                  "Female"))
## Column: SMN2 copy number
baseline_table$SMN2_CN <- 
  factor(baseline_table$SMN2_CN, levels=c(2,3,4,5),
         labels=c("2xSMN2", 
                  "3xSMN2",
                  "4xSMN2",
                  "5xSMN2"))
## Column: extended SMA type
baseline_table$SMA_type_ext <- 
  factor(baseline_table$SMA_type_ext, levels=c("1b", "1c", "2a", "2b", "3a", "3b", "4"),
         labels=c("Type 1b", 
                  "Type 1c",
                  "Type 2a",
                  "Type 2b",
                  "Type 3a", 
                  "Type 3b",
                  "Type 4"))
## Column: tissue type
baseline_table$tissue <- 
  factor(baseline_table$tissue, levels=c("blood", "fib"),
         labels=c("Blood", 
                  "Fibroblasts"))

# Add labels to baseline_table (labels are used for table printing)
label(baseline_table$SMA_type_ext) <- "SMA type"
label(baseline_table$Sex) <- "Sex"
label(baseline_table$SMN2_CN) <- "SMN2 copy number"
label(baseline_table$age_at_onset_years) <- "Age at onset"
units(baseline_table$age_at_onset_years) <- "years"
label(baseline_table$age_at_sampling) <- "Age at sampling"
units(baseline_table$age_at_sampling) <- "years"

# Generate the table
table1(~ Sex + SMA_type_ext + age_at_onset_years + age_at_sampling | tissue*SMN2_CN, data = baseline_table, overall = F)


