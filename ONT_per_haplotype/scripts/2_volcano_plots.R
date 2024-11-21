# Run the following script first: 1_process_bed_file.R

# Load required libraries for: 2_volcano_plots.R
library(tidyverse)

# Set colour palette
colour_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
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

# Make sample sheet from bed_strands_merged
mdata <- bed_strands_merged %>%
  filter(SMN1_CN == 0) %>%
  mutate(tissue_factor = dplyr::case_when(
    tissue == "blood" ~ 0L,
    tissue == "fib" ~ 1L,
    TRUE ~ NA_integer_)) %>%
  mutate(downstream_env_factor = dplyr::case_when(
    downstream_env == "SMN2" ~ 0L,
    downstream_env == "SMN1" ~ 1L,
    TRUE ~ NA_integer_))

# Return unique sites
sites <- unique(mdata$POS_ID)

# Check if model works
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

summary(check_model$model)

# BLOOD
#test if sites are significantly associated with downstream env in blood
unique(mdata$tissue)
unique(mdata$tissue_factor)

# Run run_diff_meth and check model with:
## covariates: sex, age at sampling, and
## dependent variable: percentage of methylation, and 
## independent variable: downstream_env_factor
results_lm_downstream_env_factor_blood <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "downstream_env_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"),
  data = filter(mdata, tissue == "blood"),
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Plot volcano plot
volcano_downstream_env_blood <- results_lm_downstream_env_factor_blood %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN1 vs. SMN2 env.\n(blood)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-110,110)) +
  scale_y_continuous(limits = c(0,200)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none")

volcano_downstream_env_blood

# Return number of significant sites (returns 5)
length(unique(filter(results_lm_downstream_env_factor_blood, Padj_fdr < 0.01)$site))

# Save volcano_downstream_env_blood as .svg
ggsave(
  "figures/volcano_downstream_env_blood.svg",
  plot = volcano_downstream_env_blood,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Report number of haplotypes per downstream environment in blood
length(unique(filter(filter(mdata, tissue == "blood"), downstream_env == "SMN1")$SMA_ID_tissue_hap))
length(unique(filter(filter(mdata, tissue == "blood"), downstream_env == "SMN2")$SMA_ID_tissue_hap))

# Create vector with sites that are significant
significant_sites_downstream_env_factor_blood <- unique(filter(results_lm_downstream_env_factor_blood, Padj_fdr < 0.01)$site)
significant_POS_downstream_env_factor_blood <- distinct(filter(select(mdata, POS, POS_ID), POS_ID %in% significant_sites_downstream_env_factor_blood))$POS
significant_sites_downstream_env_factor_blood
## Comment: chr5_71381516 is significant in blood and not a paraphase SNP and no SNP in our data, but only few (3) haplotypes with SMN1 env


# Plot sites that are significantly associated with downstream environment in blood
blood_significant_sites <- filter(mdata, tissue == "blood") %>%
  mutate(significant = case_when(POS_ID %in% filter(results_lm_downstream_env_factor_blood, Padj_fdr < 0.01)$site ~ "Yes",
                                 !(POS_ID %in% filter(results_lm_downstream_env_factor_blood, Padj_fdr < 0.01)$site) ~ "No")) %>%
  filter(significant == "Yes") %>%
  filter(downstream_env == "SMN1" | downstream_env == "SMN2") %>%
  ggplot(aes(x = downstream_env, y = percentage, color = downstream_env)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
  geom_point(size = 2, stroke = 0, alpha = 0.5, position = position_jitterdodge()) +
  xlab("Downstream environment") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(values = c("#F0E442", "#0072B2"), name = "Downstream env.") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(POS_ID), ncol = 7)

blood_significant_sites

# Save blood_significant_sites as .svg
ggsave(
  "figures/blood_significant_sites.svg",
  plot = blood_significant_sites,
  scale = 1,
  width = 16,
  height = 4,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# FIBROBLASTS
# test if sites are significantly associated with downstream env in fibroblasts
unique(mdata$tissue)
unique(mdata$tissue_factor)

# Run run_diff_meth and check model with:
## covariates: sex, age at sampling, and
## dependent variable: percentage of methylation, and 
## independent variable: downstream_env_factor
results_lm_downstream_env_factor_fib <- purrr::map_df(
  unique(mdata$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "downstream_env_factor",
  ID = "POS_ID",
  covar = c("Sex", "age_at_sampling", "lib_size_per_SMN_copy"), 
  data = filter(mdata, tissue == "fib"),
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Plot volcano plot (fib)
volcano_downstream_env_fib <- results_lm_downstream_env_factor_fib %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN1 vs. SMN2 env.\n(fibroblasts)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-110,110)) +
  scale_y_continuous(limits = c(0,200)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_downstream_env_fib

# Print number of significant sites (fib, returns 7)
length(unique(filter(results_lm_downstream_env_factor_fib, Padj_fdr < 0.01)$site))

# Save volcano_downstream_env_fib as .svg
ggsave(
  "figures/volcano_downstream_env_fib.svg",
  plot = volcano_downstream_env_fib,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Report number of haplotypes per downstream environment in blood
length(unique(filter(filter(mdata, tissue == "fib"), downstream_env == "SMN1")$SMA_ID_tissue_hap))
length(unique(filter(filter(mdata, tissue == "fib"), downstream_env == "SMN2")$SMA_ID_tissue_hap))

# Create vector with sites that are significant
significant_sites_downstream_env_factor_fib <- unique(filter(results_lm_downstream_env_factor_fib, Padj_fdr < 0.01)$site)
significant_POS_downstream_env_factor_fib <- distinct(filter(select(mdata, POS, POS_ID), POS_ID %in% significant_sites_downstream_env_factor_fib))$POS
significant_sites_downstream_env_factor_fib
#chr5_71422606 is significant in fibroblasts but is not a paraphase SNP, but there is a SNP present in our own ONT data

# Number of significantly associated sites that overlap between blood and fibroblasts
intersect(significant_sites_downstream_env_factor_blood, significant_sites_downstream_env_factor_fib)

# Plot sites that are significantly associated with downstream environment in fibroblasts
fib_significant_sites <- filter(mdata, tissue == "fib") %>%
  mutate(significant = case_when(POS_ID %in% filter(results_lm_downstream_env_factor_fib, Padj_fdr < 0.01)$site ~ "Yes",
                                 !(POS_ID %in% filter(results_lm_downstream_env_factor_fib, Padj_fdr < 0.01)$site) ~ "No")) %>%
  filter(significant == "Yes") %>%
  filter(downstream_env == "SMN1" | downstream_env == "SMN2") %>%
  ggplot(aes(x = downstream_env, y = percentage, color = downstream_env)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
  geom_point(size = 2, stroke = 0, alpha = 0.5, position = position_jitterdodge()) +
  xlab("Downstream environment") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(values = c("#F0E442", "#0072B2"), name = "Downstream env.") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(POS_ID), ncol = 7)
fib_significant_sites

# Save fib_significant_sites as .svg
ggsave(
  "figures/fib_significant_sites.svg",
  plot = fib_significant_sites,
  scale = 1,
  width = 16,
  height = 4,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Write results in table
write_tsv(results_lm_downstream_env_factor_blood, "results_tables/1H_results_lm_downstream_env_factor_blood.txt")
write_tsv(results_lm_downstream_env_factor_fib, "results_tables/1I_results_lm_downstream_env_factor_fib.txt")

# Clean up environment
rm(results_lm_downstream_env_factor_blood)
rm(results_lm_downstream_env_factor_fib)
rm(volcano_downstream_env_blood)
rm(volcano_downstream_env_fib)
