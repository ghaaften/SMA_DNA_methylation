# Run the following script before continuing: 1_data_merge_and_filter.R

# Load the necessary libraries
library(tidyverse)
library(table1)

# Make a new table, so the table sample_info remains as it was
baseline_table <- sample_info

# Turn columns into factors with right labels
baseline_table$sex <- 
  factor(baseline_table$sex, levels=c("M","F"),
         labels=c("Male", 
                  "Female"))
baseline_table$SMN2_CN <- 
  factor(baseline_table$SMN2_CN, levels=c(2,3,4,5),
         labels=c("2xSMN2", 
                  "3xSMN2",
                  "4xSMN2",
                  "5xSMN2"))
baseline_table$NAIP_CN <- 
  factor(baseline_table$NAIP_CN, levels=c(0,1,2,3,4),
         labels=c("0xNAIP", 
                  "1xNAIP",
                  "2xNAIP",
                  "3xNAIP",
                  "4xNAIP"))

# Add labels and units where needed for making an R table
label(baseline_table$SMA_type_simpl) <- "SMA type"
label(baseline_table$sex) <- "Sex"
label(baseline_table$SMN2_CN) <- "SMN2 copy number"
label(baseline_table$NAIP_CN) <- "NAIP copy number"
label(baseline_table$age_at_onset_years) <- "Age at onset"
units(baseline_table$age_at_onset_years) <- "years"
label(baseline_table$age_at_sampling_years) <- "Age at sampling"
units(baseline_table$age_at_sampling_years) <- "years"

# Create baseline table for samples with no SMN1 and no c.859G>C mutation
baseline_table_no_SMN1_no_859 <- baseline_table %>%
  filter(!is.na(sex)) %>%
  filter(SMN1_CN == 0) %>%
  filter(mut_859_GC == "no")

# Create baseline table for samples with SMN1
baseline_table_SMN1 <- baseline_table %>%
  filter(!is.na(sex)) %>%
  filter(SMN1_CN == 1) %>%
  filter(mut_859_GC == "no")

# Create baseline table for samples with c.859G>C mutation
baseline_table_859 <- baseline_table %>%
  filter(!is.na(sex)) %>%
  filter(SMN1_CN == 0) %>%
  filter(mut_859_GC == "yes")

# Generate the baseline table (other tables as above can be visualized with the same command, other dataset)
table1(~ sex + SMA_type_simpl + age_at_onset_years + age_at_sampling_years | SMN2_CN, data = baseline_table_no_SMN1_no_859)

#clean up variables
rm(baseline_table_no_SMN1_no_859)
rm(baseline_table_SMN1)
rm(baseline_table_859)
rm(baseline_table)

