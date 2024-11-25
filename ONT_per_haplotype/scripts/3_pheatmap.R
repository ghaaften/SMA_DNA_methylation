# Run the following scripts first: 1_process_bed_file.R and 2_volcano_plots.R

# Load required libraries
library(tidyverse)
library(pheatmap)
library(grid)

# Pivot dataframe: bed_strands_merged for pheatmap
combined_data_filtered_pivoted <- bed_strands_merged %>%
  select(Anonymized_ID_tissue_hap, start, percentage) %>%
  pivot_wider(names_from = start, values_from = percentage)

# Remove haplotypes that have sites with more than 178 NAs (more than 25% of sites)
combined_data_filtered_pivoted$na_count <- apply(is.na(combined_data_filtered_pivoted), 1, sum)

# Filter for non-NAs --> returns 72 haplotypes
combined_data_filtered_pivoted <- combined_data_filtered_pivoted %>%
  filter(na_count < 178) %>%
  select(-na_count)

# Make matrix from data frame
combined_data_filtered_pivoted <- as.data.frame(combined_data_filtered_pivoted)

## Use Anonymized_ID_tissue_hap as rowname
rownames(combined_data_filtered_pivoted) <- combined_data_filtered_pivoted$Anonymized_ID_tissue_hap

## Remove Anonymized_ID_tissue_hap from table
combined_data_filtered_pivoted <- combined_data_filtered_pivoted %>%
  select(!Anonymized_ID_tissue_hap)

## Convert dataframe to matrix
combined_data_filtered_pivoted <- as.matrix(combined_data_filtered_pivoted) #convert to matrix for pheatmap

# Make annotation table (rows)
annotation <- bed_strands_merged %>%
  select(Anonymized_ID_tissue_hap, tissue, Sex, downstream_env) %>%
  distinct()

annotation <- as.data.frame(annotation)

## Assign Anonymized_ID_tissue_hap as rowname
rownames(annotation) <- annotation$Anonymized_ID_tissue_hap

## Remove Anonymized_ID_tissue_hap from dataframe
annotation <- annotation[,-1]

## Assign new column names to annotation dataframe 
colnames(annotation) <- c("Tissue", "Sex", "Downstream environment")
unique(annotation$Tissue)
unique(annotation$Sex)
unique(annotation$`Downstream environment`)


# Make annotation table (columns) for significant sites
annotation_columns <- bed_strands_merged %>%
  mutate(blood_sig = case_when(POS %in% significant_POS_downstream_env_factor_blood ~ "yes",
                                TRUE ~ "no"),
         fib_sig = case_when(POS %in% significant_POS_downstream_env_factor_fib ~ "yes",
                               TRUE ~ "no")) %>%
  select(POS, blood_sig, fib_sig) %>% 
  distinct()

annotation_columns <- as.data.frame(annotation_columns)

## Assign POS as rowname and remove from dataframe
rownames(annotation_columns) <- annotation_columns$POS
annotation_columns <- annotation_columns[,-1]

## Rename column names
colnames(annotation_columns) <- c("Significant (blood)", "Significant (fibroblasts)")

# Define colors for annotation in pheatmap
ann_colors = list("Significant (blood)" = c("yes" = "red","no" = "white"),
                  "Significant (fibroblasts)" = c("yes" = "blue","no" = "white"),
                  Sex = c("F" = "#D55E00","M" = "#56B4E9"),
                  Tissue = c("fib" = "#E69F00",
                             "blood" = "#CC79A7"),
                  "Downstream environment" = c("SMN1" = "#F0E442", "SMN2" = "#0072B2", "mixed" = "#009E73"))

# Plot heatmap
ONT_heatmap_per_haplotype <- pheatmap(combined_data_filtered_pivoted,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean", #options: euclidean, manhattan, correlation, maximum, canberra, binary
         clustering_method = "ward.D2", #options: complete, ward.D2, single, average, mcquitty, median, centroid
         annotation_col = annotation_columns,
         annotation_row = annotation,
         annotation_names_row = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_colors = ann_colors,
         fontsize = 6
)

# Prepare for saving
ONT_heatmap_per_haplotype$gtable$grobs[[1]]$gp <- gpar(lwd = 0.5)
ONT_heatmap_per_haplotype

# Save heatmap as .svg
ggsave(
  "figures/ONT_heatmap_per_haplotype.svg",
  plot = ONT_heatmap_per_haplotype,
  scale = 1,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



