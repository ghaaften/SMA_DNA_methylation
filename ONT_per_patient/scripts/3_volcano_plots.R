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

## Run differentially methylated analysis
## Filter mdata for only Concordance, only in patients with 3 or 4 SMN2 copies (blood)
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
length(unique(filter(mdata_blood, concordance_binary == "more_severe")$SMA_ID))
# number of patients less severe (blood)
length(unique(filter(mdata_blood, concordance_binary == "less_severe")$SMA_ID))

# Concordance (fibroblasts) (differential methylation analysis)
## Run differentially methylated analysis
### Filter mdata for only Concordance, only in patients with 3 or 4 SMN2 copies (blood)
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
length(unique(filter(mdata_fib, concordance_binary == "more_severe")$SMA_ID))
## number of patients less severe (fibroblasts)
length(unique(filter(mdata_fib, concordance_binary == "less_severe")$SMA_ID))

## number of CpG sites
length(unique(mdata$POS_ID))

# Write results to tsv files
write_tsv(results_lm_tissue_factor, "results_tables/1D_results_lm_tissue_factor.txt")
write_tsv(results_lm_concordance_binary_factor_blood, "results_tables/1E_results_lm_concordance_binary_factor_blood.txt")
write_tsv(results_lm_concordance_binary_factor_fib, "results_tables/1F_results_lm_concordance_binary_factor_fib.txt")

