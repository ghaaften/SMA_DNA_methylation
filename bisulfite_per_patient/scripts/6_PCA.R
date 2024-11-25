# Required scripts to run before continuing:
# 1_data_merge_and_filter.R
# 5_volcano_and_visualisations.R

# For more information on the specific PCA methods used here, 
# see: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#google_vignette

# Load libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)

# PCA plot male vs female
mdata_sex_pivoted <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(sex == "M" | sex == "F") %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  pivot_wider(names_from = POS_ID, values_from = percentage) %>%
  left_join(distinct(select(mdata_filtered_condensed_percentage_no_SMN1_859, Anonymized_ID, sex)))

## Pivot and change rownames to prepare for PCA
mdata_sex_pivoted <- as.data.frame(mdata_sex_pivoted)
rownames(mdata_sex_pivoted) <- mdata_sex_pivoted$Anonymized_ID
mdata_sex_pivoted <- mdata_sex_pivoted %>%
  select(!Anonymized_ID)

## Running PCA and visualizing
res.pca.sex <- PCA(mdata_sex_pivoted[,-58], graph = FALSE, )
PCA_sex <- fviz_pca_ind(res.pca.sex,
                        geom.ind = "point",
                        pointsize = 0.5,
                        alpha.ind = 0.5,
                        col.ind = mdata_sex_pivoted$sex, # color by groups
                        palette = c("#E69F00", "#56B4E9"),
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Sex"
                        ) +
  scale_x_continuous(limits = c(-10,15)) +
  scale_y_continuous(limits = c(-10,10)) +
  labs(title ="PCA - Sex") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")
PCA_sex

## Saving PCA as .svg 
ggsave(
  "figures/PCA_sex.svg",
  plot = PCA_sex,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# PCA plot age group
mdata_age_group_pivoted <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(age_group == "Pediatric (<18)" | age_group == "Adult (18+)") %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  pivot_wider(names_from = POS_ID, values_from = percentage) %>%
  left_join(distinct(select(mdata_filtered_condensed_percentage_no_SMN1_859, Anonymized_ID, age_group)))

## Pivot and change rownames to prepare for PCA
mdata_age_group_pivoted <- as.data.frame(mdata_age_group_pivoted)
rownames(mdata_age_group_pivoted) <- mdata_age_group_pivoted$Anonymized_ID
mdata_age_group_pivoted <- mdata_age_group_pivoted %>%
  select(!Anonymized_ID)

## Running PCA and visualizing
res.pca.age_group <- PCA(mdata_age_group_pivoted[,-58], graph = FALSE)
PCA_age_group <- fviz_pca_ind(res.pca.age_group,
                        geom.ind = "point",
                        pointsize = 0.5,
                        alpha.ind = 0.5,
                        col.ind = mdata_age_group_pivoted$age_group, # color by groups
                        palette = c("#E69F00", "#56B4E9"),
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Age group"
                        ) +
  scale_x_continuous(limits = c(-10,15)) +
  scale_y_continuous(limits = c(-10,10)) +
  labs(title ="PCA - Age group") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")
PCA_age_group

## Saving PCA_age_group as .svg
ggsave(
  "figures/PCA_age_group.svg",
  plot = PCA_age_group,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# PCA plot SMN2 CN
mdata_SMN2_CN_pivoted <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(!is.na(SMN2_CN)) %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  pivot_wider(names_from = POS_ID, values_from = percentage) %>%
  left_join(distinct(select(mdata_filtered_condensed_percentage_no_SMN1_859, Anonymized_ID, SMN2_CN)))

## Pivot dataframe and change rownames
mdata_SMN2_CN_pivoted <- as.data.frame(mdata_SMN2_CN_pivoted)
rownames(mdata_SMN2_CN_pivoted) <- mdata_SMN2_CN_pivoted$Anonymized_ID
mdata_SMN2_CN_pivoted <- mdata_SMN2_CN_pivoted %>%
  select(!Anonymized_ID)

## Run PCA and visualize
res.pca.SMN2_CN <- PCA(mdata_SMN2_CN_pivoted[,-58], graph = FALSE)
PCA_SMN2_CN <- fviz_pca_ind(res.pca.SMN2_CN,
                              geom.ind = "point",
                              pointsize = 0.5,
                              alpha.ind = 0.5,
                              col.ind = as.factor(mdata_SMN2_CN_pivoted$SMN2_CN), # color by groups
                              palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                              addEllipses = TRUE, # Concentration ellipses
                              legend.title = "SMN2 CN"
                            ) +
  scale_x_continuous(limits = c(-10,15)) +
  scale_y_continuous(limits = c(-11,9)) +
  labs(title ="PCA - SMN2 CN") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")
PCA_SMN2_CN

## Save PCA_SMN2_CN as .svg
ggsave(
  "figures/PCA_SMN2_CN.svg",
  plot = PCA_SMN2_CN,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# PCA plot SMA_type
mdata_SMA_type_pivoted <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(!is.na(SMA_type_simpl)) %>%
  filter(SMA_type_simpl != "Presymptomatic") %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  pivot_wider(names_from = POS_ID, values_from = percentage) %>%
  left_join(distinct(select(mdata_filtered_condensed_percentage_no_SMN1_859, Anonymized_ID, SMA_type_simpl)))

## Pivot dataframe and change rownames
mdata_SMA_type_pivoted <- as.data.frame(mdata_SMA_type_pivoted)
rownames(mdata_SMA_type_pivoted) <- mdata_SMA_type_pivoted$Anonymized_ID
mdata_SMA_type_pivoted <- mdata_SMA_type_pivoted %>%
  select(!Anonymized_ID)

## Run PCA and visualize
res.pca.SMA_type <- PCA(mdata_SMA_type_pivoted[,-58], graph = FALSE)
PCA_SMA_type <- fviz_pca_ind(res.pca.SMA_type,
                            geom.ind = "point",
                            pointsize = 0.5,
                            alpha.ind = 0.5,
                            col.ind = mdata_SMA_type_pivoted$SMA_type_simpl, # color by groups
                            palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"),
                            addEllipses = TRUE, # Concentration ellipses
                            legend.title = "SMA type"
                            ) +
  scale_x_continuous(limits = c(-10,15)) +
  scale_y_continuous(limits = c(-11,9)) +
  labs(title ="PCA - SMA type") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")
PCA_SMA_type

## Save PCA_SMA_type as .svg
ggsave(
  "figures/PCA_SMA_type.svg",
  plot = PCA_SMA_type,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Treatment
## PCA plot treatment groups HFMSE
## keep only 3 and 4 copy groups because 2 and 5 copy groups are very small
mdata_filtered_condensed_percentage_hfmse <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4)
mdata_filtered_condensed_percentage_hfmse_pivoted <- mdata_filtered_condensed_percentage_hfmse %>%
  filter(!is.na(dHFMSE)) %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  pivot_wider(names_from = POS_ID, values_from = percentage) %>%
  left_join(distinct(select(mdata_filtered_condensed_percentage_hfmse, Anonymized_ID, HFMSE_response_group)))

mdata_filtered_condensed_percentage_hfmse_pivoted <- as.data.frame(mdata_filtered_condensed_percentage_hfmse_pivoted)
rownames(mdata_filtered_condensed_percentage_hfmse_pivoted) <- mdata_filtered_condensed_percentage_hfmse_pivoted$Anonymized_ID
mdata_filtered_condensed_percentage_hfmse_pivoted <- mdata_filtered_condensed_percentage_hfmse_pivoted %>%
  select(!Anonymized_ID)

## Run PCA
res.pca.dHFMSE <- PCA(mdata_filtered_condensed_percentage_hfmse_pivoted[,-58], graph = FALSE)
PCA_dHFMSE <- fviz_pca_ind(res.pca.dHFMSE,
                           geom.ind = "point",
                           pointsize = 0.5,
                           alpha.ind = 0.5,
                           col.ind = mdata_filtered_condensed_percentage_hfmse_pivoted$HFMSE_response_group, # color by groups
                           palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"),
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "dHFMSE"
                           ) +
  labs(title ="PCA - dHFMSE") +
  scale_x_continuous(limits = c(-9,22)) +
  scale_y_continuous(limits = c(-12,9)) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8), legend.position = "none")
PCA_dHFMSE

## Save PCA_dHFMSE as .svg
ggsave(
  "figures/PCA_dHFMSE.svg",
  plot = PCA_dHFMSE,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

