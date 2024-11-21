# Run the following script before continuing: 1_data_merge_and_filter.R

# Load required libraries for "5_PCA_plots.R"
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)


# PCA plot standards vs SMA
## Pivot bed_strands_merged
PCA_pivoted <- bed_strands_merged %>%
  select(SMA_ID_tissue, POS, percentage) %>%
  pivot_wider(names_from = POS, values_from = percentage) %>%
  left_join(distinct(select(bed_strands_merged, SMA_ID_tissue, tissue))) 

PCA_pivoted <- as.data.frame(PCA_pivoted)

## Use SMA_ID_tissue as rownames
rownames(PCA_pivoted) <- PCA_pivoted$SMA_ID_tissue

## Remove SMA_ID_tissue column
PCA_pivoted <- PCA_pivoted %>%
  select(!SMA_ID_tissue)

## Make PCA analysis (tissue)
res.pca <- PCA(PCA_pivoted[,-711], graph = FALSE)

## Plot res.pca
PCA_tissues <- fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointsize = 0.5,
             alpha.ind = 0.5,
             linewidth = 0.25,
             col.ind = PCA_pivoted$tissue,
             palette = c("#CC79A7", "#E69F00"), #
             addEllipses = TRUE, 
             legend.title = "Tissue type") +
  labs(title ="PCA - Tissue type") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")

PCA_tissues

## Save PCA_tissues as .svg
ggsave(
  "figures/PCA_tissues.svg",
  plot = PCA_tissues,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

