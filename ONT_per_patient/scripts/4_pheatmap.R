# Run the following scripts before continuing: 1_data_merge_and_filter.R and 3_volcano_plots.R

# Load required libraries for "4_pheatmap.R"
library(tidyverse)
library(pheatmap)
library(grid)

# Pivot dataframe bed_strands_merged for pheatmap
combined_data_filtered_pivoted <- bed_strands_merged %>%
  select(SMA_ID_tissue, start, percentage) %>%
  pivot_wider(names_from = start, values_from = percentage)

# Return as data frame
combined_data_filtered_pivoted <- as.data.frame(combined_data_filtered_pivoted)

# Use SMA_ID_tissue as rownames
rownames(combined_data_filtered_pivoted) <- combined_data_filtered_pivoted$SMA_ID_tissue

# Remove SMA_ID_tissue column, as that column is now used as rowname
combined_data_filtered_pivoted <- combined_data_filtered_pivoted %>%
  select(!SMA_ID_tissue)

# Convert df into matrix
combined_data_filtered_pivoted <- as.matrix(combined_data_filtered_pivoted)

# Make annotation table for pheatmap (rows)
annotation <- bed_strands_merged %>%
  select(SMA_ID_tissue, tissue, Sex, concordance_binary) %>%
  distinct()

annotation <- as.data.frame(annotation)

# Set SMA_ID_tissue as rownames
rownames(annotation) <- annotation$SMA_ID_tissue

# Remove SMA_ID_tissue from annotation df
annotation <- annotation[,-1]

# Rename columns of annotation df
colnames(annotation) <- c("Tissue", "Sex", "Concordance")

# Check unique values of these columns
unique(annotation$Tissue)
unique(annotation$Sex)
unique(annotation$Concordance)

# Make annotation table for pheatmap (columns)
annotation_columns <- bed_strands_merged %>%
  mutate(tissue_sig = case_when(POS %in% significant_POS_tissue ~ "yes",
                                TRUE ~ "no")) %>%
  select(POS, Bisulfite_amplicon, tissue_sig) %>% 
  distinct()

annotation_columns <- as.data.frame(annotation_columns)

# Set rownames as POS (positions)
rownames(annotation_columns) <- annotation_columns$POS

# Remove SMA_ID_tissue from annotation_columns df
annotation_columns <- annotation_columns[,-1]

# Rename columns of annotation_columns df
colnames(annotation_columns) <- c("Bisulfite_amplicon", "Tissue_significant")



# Define colors for annotation in pheatmap
ann_colors = list(Bisulfite_amplicon = c("yes" = "blue",
                                         "no" = "white"),
                  Tissue_significant = c("yes" = "red",
                                         "no" = "white"),
                  Tissue = c("fib" = "#E69F00",
                          "blood" = "#CC79A7"),
                  Sex = c("F" = "#D55E00",
                          "M" = "#56B4E9"),
                  Concordance = c("more_severe" = "#F0E442",
                                  "less_severe" = "#009E73"))

# Make heatmap with pheatmap (per patient)
## All SMN positions and all SMA samples
ONT_heatmap_per_patient <- pheatmap(combined_data_filtered_pivoted,
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
         fontsize = 6, 
)

# Adjust line width of dendrogram 
ONT_heatmap_per_patient$gtable$grobs[[1]]$gp <- gpar(lwd = 0.5)
ONT_heatmap_per_patient

# Save figure ONT_heatmap_per_patient as .svg
ggsave(
  "figures/ONT_heatmap_per_patient.svg",
  plot = ONT_heatmap_per_patient,
  scale = 1,
  width = 16,
  height = 8,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

