# Run the following script before continuing: 1_data_merge_and_filter.R

# Load required libraries for "3_volcano_plots.R"
library(tidyverse)

# Define colour palette
colour_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#000000",
                    "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

# Define function to run differential methylation analysis
run_diff_meth <- function (site, 
                           dependent_variable, 
                           independent_variable, 
                           data, 
                           ID = "POS_ID", 
                           covar = NULL, 
                           returnModel = FALSE,
                           method = c("lm", "glm")
) {
  
  # specify method 
  method <- match.arg(method)
  
  # subset site
  dat <- data[data[[ID]] == site,]
  
  # nr of samples with non-missing dependent,independent and covar
  nr_samples <- sum(complete.cases(dat[,c(dependent_variable, independent_variable, covar), drop = FALSE]))
  
  # Define formula for model to run (including or excluding covariates)
  if (!is.null(covar)) {
    # Including covariates
    formula <- as.formula(sprintf("%s ~ %s + %s", dependent_variable, independent_variable, paste(covar, collapse = "+")))
    formula_without_independent_variable <- as.formula(sprintf("%s ~ %s", dependent_variable, paste(covar, collapse = "+")))
  } else {
    # Excluding covariates
    formula <- as.formula(sprintf("%s ~ %s", dependent_variable, independent_variable))
    formula_without_independent_variable <- as.formula(sprintf("%s ~ 1", dependent_variable))
  }
  
  # Make model based on defined formula
  if (method == "lm") {
    # Linear model
    model <- lm(formula, data = dat)
    model_without_independent_variable <- lm(formula_without_independent_variable, data = dat)
  } else if (method == "glm") {
    # Generalized linear model
    model <- glm(formula, data = dat, family = "binomial")
    model_without_independent_variable <- glm(formula_without_independent_variable, data = dat, family = "binomial")
  }
  
  # Create test statistics dataframe
  stats <- tibble(
    site = site,
    independent_variable = independent_variable,
    dependent_variable = dependent_variable,
    N = nr_samples,
    b = summary(model)$coef[independent_variable, 1],
    se = summary(model)$coef[independent_variable, 2], 
    b_95CI_lower = confint(model)[independent_variable, 1],
    b_95CI_upper = confint(model)[independent_variable, 2],
    t = summary(model)$coef[independent_variable, 3], 
    P = summary(model)$coef[independent_variable, 4],
    R2AB = summary(model)$r.squared, # R squared of model with independent variable and covariates
    R2A = summary(model_without_independent_variable)$r.squared # R squared of model with only covariates
  ) %>%
    mutate(f2 = (R2AB - R2A) / (1 - R2AB)) #Cohen's f2, local effect size
  
  # Return stats dataframe and model (if returnModel = TRUE)
  if (returnModel) {
    list(
      stats = stats,
      model = if(returnModel) model else NULL
    )
  } else {
    stats
  }
}


# Sample sheet bed_strands_merged (script "1_data_merge_and_filter.R")
mdata <- bed_strands_merged %>%
  filter(SMN1_CN == 0) %>%
  mutate(tissue_factor = dplyr::case_when(
    tissue == "blood" ~ 0L,
    tissue == "fib" ~ 1L,
    TRUE ~ NA_integer_)) %>%
  mutate(concordance_binary_factor = dplyr::case_when(
    # Add binary factor for disease concordance
    concordance_binary == "less_severe" ~ 0L,
    concordance_binary == "more_severe" ~ 1L,
    TRUE ~ NA_integer_))

## Filter for unique sites
sites <- unique(mdata$POS_ID)

# Run run_diff_meth and check model with:
## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: tissue type
check_model <- run_diff_meth(
  site = "chr5_71382829",
  dependent_variable = "percentage",
  independent_variable = "tissue_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata,
  returnModel = TRUE,
  method = "lm"
)

## Show summary
summary(check_model$model)

# Tissue (differential methylation analysis)

## Test if sites are significantly associated with tissue
unique(mdata$tissue)
unique(mdata$tissue_factor)

## Run run_diff_meth function with 'lm' with:
## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: tissue type
results_lm_tissue_factor <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "tissue_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Visualize differential methylation with volcano plot
volcano_tissue <- results_lm_tissue_factor %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Tissue type") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-60,60)) +
  scale_y_continuous(limits = c(0,15)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none")

volcano_tissue

# Check number of differentially methylated sites between blood and fibroblasts (returns 58)
length(unique(filter(results_lm_tissue_factor, Padj_fdr < 0.01)$site))

ggsave(
  "figures/volcano_tissue.svg",
  plot = volcano_tissue,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Concordance (blood) (differential methylation analysis)

## Make a vector with sites that are significant
significant_sites_tissue <- unique(filter(results_lm_tissue_factor, Padj_fdr < 0.01)$site)
significant_POS_tissue <- distinct(filter(select(bed_strands_merged, POS, POS_ID), POS_ID %in% significant_sites_tissue))$POS

## Run differential methylation analysis
## Filter mdata for only patients with 3 or 4 SMN2 copies (blood)
mdata_blood <- mdata %>%
  filter(tissue == "blood") %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4)

## Run run_diff_meth function with 'lm' with:
## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: concordance_binary_factor
results_lm_concordance_binary_factor_blood <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "concordance_binary_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"), 
  data = mdata_blood,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))


## Visualizing volcano plot
volcano_concordance_binary_factor_blood <- results_lm_concordance_binary_factor_blood %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Concordance (blood)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-60,60)) +
  scale_y_continuous(limits = c(0,15)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none")

volcano_concordance_binary_factor_blood

## Return length of significant sites (returns: 0)
length(unique(filter(results_lm_concordance_binary_factor_blood, Padj_fdr < 0.01)$site))

## Save volcano_concordance_binary_factor_blood figure as .svg
ggsave(
  "figures/volcano_concordance_binary_factor_blood.svg",
  plot = volcano_concordance_binary_factor_blood,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# number of patients more severe (blood)
length(unique(filter(mdata_blood, concordance_binary == "more_severe")$Anonymized_ID))
# number of patients less severe (blood)
length(unique(filter(mdata_blood, concordance_binary == "less_severe")$Anonymized_ID))

# Concordance (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 3 or 4 SMN2 copies (blood)
mdata_fib <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: concordance_binary_factor
results_lm_concordance_binary_factor_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "concordance_binary_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_concordance_binary_factor_fib <- results_lm_concordance_binary_factor_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Concordance (fibroblasts)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-60,60)) +
  scale_y_continuous(limits = c(0,15)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
        )

volcano_concordance_binary_factor_fib

# Return number of significant sites
length(unique(filter(results_lm_concordance_binary_factor_fib, Padj_fdr < 0.01)$site))

## Save volcano_concordance_binary_factor_fib figure as .svg
ggsave(
  "figures/volcano_concordance_binary_factor_fib.svg",
  plot = volcano_concordance_binary_factor_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients more severe (fibroblasts)
length(unique(filter(mdata_fib, concordance_binary == "more_severe")$Anonymized_ID))
## number of patients less severe (fibroblasts)
length(unique(filter(mdata_fib, concordance_binary == "less_severe")$Anonymized_ID))



# Plot SMN2-FL RNA expression with age
RNA_SMN2_age_graph <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(SMN2_ave)) %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, SMN2_ave) %>%
  distinct() %>%
  ggplot(aes(x = age_at_sampling, y = SMN2_ave, colour = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_point(size = 1, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  xlab("Age at sampling (years)") +
  ylab("SMN2-FL expression (a.u.)") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(SMN2_CN), ncol = 2)

RNA_SMN2_age_graph

# Save RNA_SMN2_age_graph as .svg
ggsave(
  "figures/RNA_SMN2_age.svg",
  plot = RNA_SMN2_age_graph,
  scale = 1,
  width = 8,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Statistics: is SMN2-FL RNA expression in fibroblasts associated with age (per copy number group)
RNA_SMN2_age_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(SMN2_ave)) %>%
  filter(SMN2_CN == 3 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, SMN2_ave) %>%
  distinct()
model_SMN2_age_CN3 <- lm(SMN2_ave ~ age_at_sampling, data = RNA_SMN2_age_CN3)
summary(model_SMN2_age_CN3)
rm(RNA_SMN2_age_CN3)
rm(model_SMN2_age_CN3)

RNA_SMN2_age_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(SMN2_ave)) %>%
  filter(SMN2_CN == 4 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, SMN2_ave) %>%
  distinct()
model_SMN2_age_CN4 <- lm(SMN2_ave ~ age_at_sampling, data = RNA_SMN2_age_CN4)
summary(model_SMN2_age_CN4)
rm(RNA_SMN2_age_CN4)
rm(model_SMN2_age_CN4)

# Plot D7 RNA expression with age
RNA_D7_age_graph <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(D7_ave)) %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, D7_ave) %>%
  distinct() %>%
  ggplot(aes(x = age_at_sampling, y = D7_ave, colour = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_point(size = 1, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  xlab("Age at sampling (years)") +
  ylab("SMN2∆7 RNA expression (a.u.)") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  theme_classic(base_size = 7.5) + 
  #theme(legend.position = "none") +
  facet_wrap(vars(SMN2_CN), ncol = 2)

RNA_D7_age_graph

# Save RNA_D7_age_graph as .svg
ggsave(
  "figures/RNA_D7_age.svg",
  plot = RNA_D7_age_graph,
  scale = 1,
  width = 8,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Statistics: is SMN2-D7 RNA expression in fibroblasts associated with age (per copy number group)
RNA_D7_age_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(D7_ave)) %>%
  filter(SMN2_CN == 3 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, D7_ave) %>%
  distinct()
model_D7_age_CN3 <- lm(D7_ave ~ age_at_sampling, data = RNA_D7_age_CN3)
summary(model_D7_age_CN3)
rm(RNA_D7_age_CN3)
rm(model_D7_age_CN3)

RNA_D7_age_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(D7_ave)) %>%
  filter(SMN2_CN == 4 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, D7_ave) %>%
  distinct()
model_D7_age_CN4 <- lm(D7_ave ~ age_at_sampling, data = RNA_D7_age_CN4)
summary(model_D7_age_CN4)
rm(RNA_D7_age_CN4)
rm(model_D7_age_CN4)

# Plot AS1 RNA expression with age
RNA_AS1_age_graph <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(AS1_ave)) %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, AS1_ave) %>%
  distinct() %>%
  ggplot(aes(x = age_at_sampling, y = AS1_ave, colour = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_point(size = 1, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  xlab("Age at sampling (years)") +
  ylab("SMN-AS1 RNA expression (a.u.)") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  theme_classic(base_size = 7.5) + 
  #theme(legend.position = "none") +
  facet_wrap(vars(SMN2_CN), ncol = 2)

RNA_AS1_age_graph

# Save RNA_AS1_age_graph as .svg
ggsave(
  "figures/RNA_AS1_age.svg",
  plot = RNA_AS1_age_graph,
  scale = 1,
  width = 8,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Statistics: is SMN-AS1 RNA expression in fibroblasts associated with age (per copy number group)
RNA_AS1_age_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(AS1_ave)) %>%
  filter(SMN2_CN == 3 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, AS1_ave) %>%
  distinct()
model_AS1_age_CN3 <- lm(AS1_ave ~ age_at_sampling, data = RNA_AS1_age_CN3)
summary(model_AS1_age_CN3)
rm(RNA_AS1_age_CN3)
rm(model_AS1_age_CN3)

RNA_AS1_age_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(!is.na(AS1_ave)) %>%
  filter(SMN2_CN == 4 ) %>% 
  select(Anonymized_ID, SMN2_CN, age_at_sampling, AS1_ave) %>%
  distinct()
model_AS1_age_CN4 <- lm(AS1_ave ~ age_at_sampling, data = RNA_AS1_age_CN4)
summary(model_AS1_age_CN4)
rm(RNA_AS1_age_CN4)
rm(model_AS1_age_CN4)

#####

# SMN2 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 3 copies (fibroblasts)
mdata_fib_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 3)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: SMN2_ave
results_lm_SMN2_RNA_CN3_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMN2_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN3,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_SMN2_RNA_CN3_fib <- results_lm_SMN2_RNA_CN3_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2-FL RNA expression\n(3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-85,85)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_SMN2_RNA_CN3_fib

# Return number of significant sites; returns 0
length(unique(filter(results_lm_SMN2_RNA_CN3_fib, Padj_fdr < 0.01)$site))

## Save volcano_SMN2_RNA_CN3_fib figure as .svg
ggsave(
  "figures/volcano_SMN2_RNA_CN3_fib.svg",
  plot = volcano_SMN2_RNA_CN3_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 3xSMN2 in fib group, that have SMN2 RNA data available; returns 10
length(unique(filter(mdata_fib_CN3, !is.na(SMN2_ave))$Anonymized_ID))




# SMN2 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 4 copies (fibroblasts)
mdata_fib_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 4)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: SMN2_ave
results_lm_SMN2_RNA_CN4_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMN2_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN4,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_SMN2_RNA_CN4_fib <- results_lm_SMN2_RNA_CN4_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2-FL RNA expression\n(4xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-85,85)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_SMN2_RNA_CN4_fib

# Return number of significant sites; returns 0
length(unique(filter(results_lm_SMN2_RNA_CN4_fib, Padj_fdr < 0.01)$site))

## Save volcano_SMN2_RNA_CN4_fib figure as .svg
ggsave(
  "figures/volcano_SMN2_RNA_CN4_fib.svg",
  plot = volcano_SMN2_RNA_CN4_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 4xSMN2 in fib group that have SMN2 RNA data available; returns 8
length(unique(filter(mdata_fib_CN4, !is.na(SMN2_ave))$Anonymized_ID))



# D7 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 3 copies (fibroblasts)
mdata_fib_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 3)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: D7_ave
results_lm_D7_RNA_CN3_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "D7_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN3,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_D7_RNA_CN3_fib <- results_lm_D7_RNA_CN3_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2∆7 RNA expression\n(3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-85,85)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_D7_RNA_CN3_fib

# Return number of significant sites; returns 0
length(unique(filter(results_lm_D7_RNA_CN3_fib, Padj_fdr < 0.01)$site))

## Save volcano_D7_RNA_CN3_fib figure as .svg
ggsave(
  "figures/volcano_D7_RNA_CN3_fib.svg",
  plot = volcano_D7_RNA_CN3_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 3xSMN2 in fib group that have D7 RNA data available; returns 10
length(unique(filter(mdata_fib_CN3, !is.na(D7_ave))$Anonymized_ID))


# D7 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 4 copies (fibroblasts)
mdata_fib_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 4)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: D7_ave
results_lm_D7_RNA_CN4_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "D7_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN4,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_D7_RNA_CN4_fib <- results_lm_D7_RNA_CN4_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2∆7 RNA expression\n(4xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-85,85)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_D7_RNA_CN4_fib

# Return number of significant sites; returns 1
length(unique(filter(results_lm_D7_RNA_CN4_fib, Padj_fdr < 0.01)$site))

## Save volcano_D7_RNA_CN4_fib figure as .svg
ggsave(
  "figures/volcano_D7_RNA_CN4_fib.svg",
  plot = volcano_D7_RNA_CN4_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 4xSMN2 in fib group that have D7 RNA data available; returns 8
length(unique(filter(mdata_fib_CN4, !is.na(D7_ave))$Anonymized_ID))

# plot significant site
plot_D7_expression <- mdata_fib %>%
  filter(POS_ID %in% filter(results_lm_D7_RNA_CN4_fib, Padj_fdr < 0.01)$site) %>%
  filter(!is.na(D7_ave)) %>%
  ggplot(aes(x = D7_ave, y = percentage, color = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_smooth(method = "lm", linewidth = 0.5) + #color = "black",
  geom_point(size = 1, stroke = 0, alpha = 0.8) +
  scale_y_continuous(limits = c(0,100)) +
  ggtitle("CpG site chr5:71417351") +
  xlab("SMN2∆7 RNA expression (a.u.)") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "SMN2 copy number") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "SMN2 copy number") +
  theme_classic(base_size = 7.5)
plot_D7_expression

# Save sites plot as .svg
ggsave(
  "figures/plot_D7_expression.svg",
  plot = plot_D7_expression,
  scale = 1,
  width = 7,
  height = 5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



# AS1 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 3 copies (fibroblasts)
mdata_fib_CN3 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 3)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: AS1_ave
results_lm_AS1_RNA_CN3_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "AS1_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN3,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_AS1_RNA_CN3_fib <- results_lm_AS1_RNA_CN3_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN-AS1 RNA expression\n(3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-6000,6000)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_AS1_RNA_CN3_fib

# Return number of significant sites; returns 0
length(unique(filter(results_lm_AS1_RNA_CN3_fib, Padj_fdr < 0.01)$site))

## Save volcano_AS1_RNA_CN3_fib figure as .svg
ggsave(
  "figures/volcano_AS1_RNA_CN3_fib.svg",
  plot = volcano_AS1_RNA_CN3_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 3xSMN2 in fib group that have AS1 RNA data available; returns 9
length(unique(filter(mdata_fib_CN3, !is.na(AS1_ave))$Anonymized_ID))



# AS1 RNA expression (fibroblasts) (differential methylation analysis)
## Run differential methylation analysis
### Filter mdata for only patients with 4 copies (fibroblasts)
mdata_fib_CN4 <- mdata %>%
  filter(tissue == "fib") %>%
  filter(SMN2_CN == 4)

## covariates: sex, age at sampling and library size, 
## dependent variable: percentage of methylation, and 
## independent variable: AS1_ave
results_lm_AS1_RNA_CN4_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "AS1_ave",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = mdata_fib_CN4,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

## Visualizing volcano plot
volcano_AS1_RNA_CN4_fib <- results_lm_AS1_RNA_CN4_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN-AS1 RNA expression\n(4xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-6000,6000)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_AS1_RNA_CN4_fib

# Return number of significant sites; returns 0
length(unique(filter(results_lm_AS1_RNA_CN4_fib, Padj_fdr < 0.01)$site))

## Save volcano_AS1_RNA_CN4_fib figure as .svg
ggsave(
  "figures/volcano_AS1_RNA_CN4_fib.svg",
  plot = volcano_AS1_RNA_CN4_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

## number of patients with 4xSMN2 in fib group that have AS1 RNA data available; returns 7
length(unique(filter(mdata_fib_CN4, !is.na(AS1_ave))$Anonymized_ID))


# Summary statistics, only fibroblast
summary_stats_fib <- mdata %>%
  filter(tissue == "fib") %>%
  group_by(POS) %>%
  summarise(mean_percentage = mean(percentage),
            stdev_percentage = sd(percentage),
            median_percentage = median(percentage),
            IQR_percentage = IQR(percentage),
            min_percentage = min(percentage),
            max_percentage = max(percentage)
  )

# Summary statistics, only blood
summary_stats_blood <- mdata %>%
  filter(tissue == "blood") %>%
  group_by(POS) %>%
  summarise(mean_percentage = mean(percentage),
            stdev_percentage = sd(percentage),
            median_percentage = median(percentage),
            IQR_percentage = IQR(percentage),
            min_percentage = min(percentage),
            max_percentage = max(percentage)
  )

# Plot mean methylation percentage for fibroblasts
mean_percentage_fib <- summary_stats_fib %>%
  filter(between(POS, 71380000, 71412000)) %>%
  ggplot(aes(x = POS, y = mean_percentage)) +
  geom_smooth(method = "loess", span = 0.05, colour = "#E69F00") +
  scale_x_continuous(limits = c(71380000, 71412000)) +
  scale_y_continuous(limits = c(-10,110), breaks = c(0,50,100)) +
  xlab("Coordinate") +
  ylab("Mean %") +
  theme_classic()
mean_percentage_fib

ggsave(
  "figures/mean_percentage_fib.svg",
  plot = mean_percentage_fib,
  scale = 1,
  width = 16,
  height = 2,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Plot mean methylation percentage for blood
mean_percentage_blood <- summary_stats_blood %>%
  filter(between(POS, 71380000, 71412000)) %>%
  ggplot(aes(x = POS, y = mean_percentage)) +
  geom_smooth(method = "loess", span = 0.05, colour = "#CC79A7") +
  scale_x_continuous(limits = c(71380000, 71412000)) +
  scale_y_continuous(limits = c(-10,110), breaks = c(0,50,100)) +
  xlab("Coordinate") +
  ylab("Mean %") +
  theme_classic()
mean_percentage_blood

ggsave(
  "figures/mean_percentage_blood.svg",
  plot = mean_percentage_blood,
  scale = 1,
  width = 16,
  height = 2,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)




## number of CpG sites
length(unique(mdata$POS_ID))

# Write summary statistics to tsv file
write_tsv(summary_stats_blood, "results_tables/summary_stats_blood.txt")
write_tsv(summary_stats_fib, "results_tables/summary_stats_fib.txt")

# Write results to tsv files
write_tsv(results_lm_tissue_factor, "results_tables/1D_results_lm_tissue_factor.txt")
write_tsv(results_lm_concordance_binary_factor_blood, "results_tables/1E_results_lm_concordance_binary_factor_blood.txt")
write_tsv(results_lm_concordance_binary_factor_fib, "results_tables/1F_results_lm_concordance_binary_factor_fib.txt")
write_tsv(results_lm_SMN2_RNA_CN3_fib, "results_tables/S2B_results_lm_SMN2_RNA_CN3_fib.txt")
write_tsv(results_lm_SMN2_RNA_CN4_fib, "results_tables/S2C_results_lm_SMN2_RNA_CN4_fib.txt")
write_tsv(results_lm_D7_RNA_CN3_fib, "results_tables/S2E_results_lm_D7_RNA_CN3_fib.txt")
write_tsv(results_lm_D7_RNA_CN4_fib, "results_tables/S2F_results_lm_D7_RNA_CN4_fib.txt")
write_tsv(results_lm_AS1_RNA_CN3_fib, "results_tables/S2I_results_lm_AS1_RNA_CN3_fib.txt")
write_tsv(results_lm_AS1_RNA_CN4_fib, "results_tables/S2J_results_lm_AS1_RNA_CN4_fib.txt")




