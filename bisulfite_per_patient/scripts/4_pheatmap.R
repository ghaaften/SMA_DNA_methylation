# Run the following script before continuing: 1_data_merge_and_filter.R

# Load required libraries to run 4_pheatmap.R
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(grid)

# Pivot the combined_filtered dataframe to make it suitable for the pheatmap tool.
combined_data_filtered_pivoted <- combined_data_filtered %>%
  select(SMA_ID, start, percentage) %>%
  pivot_wider(names_from = start, values_from = percentage)

combined_data_filtered_pivoted <- as.data.frame(combined_data_filtered_pivoted)

# Set SMA_ID as rowname and remove the column from dataframe
rownames(combined_data_filtered_pivoted) <- combined_data_filtered_pivoted$SMA_ID
combined_data_filtered_pivoted <- combined_data_filtered_pivoted %>%
  select(!SMA_ID)
combined_data_filtered_pivoted <- as.matrix(combined_data_filtered_pivoted) #convert to matrix for pheatmap

# Create a dataframe that can be used for annotation for rownames in pheatmap
annotation <- sample_info %>%
  select(SMA_ID, sex, SMN2_CN, SMA_type_simpl, age_at_sampling_years) #  age_at_sampling_years, age_at_onset_years
annotation <- as.data.frame(annotation)
rownames(annotation) <- annotation$SMA_ID
annotation <- annotation[,-1]
colnames(annotation) <- c("Sex", "SMN2 copy number", "SMA type", "Age at sampling (years)")

# Define colors for row annotation in pheatmap
ann_colors = list("Sex" = c("F" = "#D55E00",
                           "M" = "#56B4E9"),
                  "SMN2 copy number" = c("2" = "#B8A92F",
                               "3" = "#E2C842",
                               "4" = "#F6EB7D",
                               "5" = "#FFF9CC"),
                  "SMA type" = c("Type 1" = "#5d133c",
                                "Type 2" = "#824163",
                                "Type 3" = "#a86e8c",
                                "Type 4" = "#cf9bb7",
                                "Presymptomatic" = "#f6cbe3"))

# Visualize the heatmap with pheatmap tool
bisulfite_heatmap <- pheatmap(combined_data_filtered_pivoted,
         border_color = NA,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean", #options: euclidean, manhattan, correlation, maximum, canberra, binary
         clustering_method = "ward.D2", #options: complete, ward.D2, single, average, mcquitty, median, centroid
         annotation_row = annotation,
         annotation_names_row = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_colors = ann_colors, #give custom colours to annotations
         fontsize = 6,
         gaps_col = c(61, 87, 100, 112, 121) #separate the different regions with whitespace: promoter, lncRNA, intron1, intron2a, 3'UTR
         )
bisulfite_heatmap$gtable$grobs[[1]]$gp <- gpar(lwd = 0.5)
bisulfite_heatmap

# Save heatmap as .svg
ggsave(
  "figures/bisulfite_heatmap.svg",
  plot = bisulfite_heatmap,
  scale = 1,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

#clean up variables
rm(annotation)
rm(ann_colors)
rm(bisulfite_heatmap)
rm(combined_data_filtered_pivoted)


