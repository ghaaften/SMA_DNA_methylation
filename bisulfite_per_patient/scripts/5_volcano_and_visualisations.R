# Run the following script before continuing: 1_data_merge_and_filter.R

# Load required libraries
library(tidyverse)
library(corrplot)
library(svglite)

# Make function to run differential methylation analysis
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

# Define colour palette
colour_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#000000",
                    "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

# Make a new variable mdata for the methylation data, so the original combined_data_filtered remains unchanged
mdata <- combined_data_filtered
sites <- unique(mdata$POS_ID)

# Calculate standard deviation (SD) of methylation percentage per site
mdata_SD_per_site <- mdata %>%
  group_by(POS_ID) %>%
  summarise(stdev_percentage_per_site = sd(percentage))

# Add extra variable in mdata that indicates whether methylation SD is lower or higher than 5
mdata <- mdata %>%
  left_join(mdata_SD_per_site) %>%
  mutate(percentage_SD_lower_5 = case_when(stdev_percentage_per_site < 5 ~ "yes",
                                         stdev_percentage_per_site >= 5 ~ "no"))
## how many sites have percentage SD < 5?
length(unique(filter(mdata, percentage_SD_lower_5 == "yes")$POS_ID)) #57 sites low SD
length(unique(filter(mdata, percentage_SD_lower_5 == "no")$POS_ID)) #64 sites high SD -> use this as multiple testing factor

## filter out sites with percentage SD lower than 5
mdata_filtered <- mdata %>%
  filter(percentage_SD_lower_5 == "no")

# Correlation plot
## Restructure the data and make a correlation plot, to determine the correlation between all CpG sites
mdata_filtered_sd_spread <- mdata_filtered %>%
  select(Anonymized_ID, POS_ID, percentage) %>%
  spread(key = POS_ID, value = percentage)

## Remove NAs
mdata_filtered_sd_spread_na_omit <- na.omit(mdata_filtered_sd_spread)
rownames(mdata_filtered_sd_spread_na_omit) <- mdata_filtered_sd_spread_na_omit[,1]
mdata_filtered_sd_spread_na_omit <- mdata_filtered_sd_spread_na_omit[,-1]
mdata_filtered_sd_spread_na_omit %>% cor(method = "spearman") -> sites_corplot
correlation_plot <- corrplot(sites_corplot, method="circle", type = 'lower', tl.col = 'black', tl.cex = 0.5,
         cl.ratio = 0.2, tl.srt = 45)

correlation_plot

## Check the number of sites 
sites_corplot_0.9 <- subset(as.data.frame.table(sites_corplot), Freq > 0.9)
sites_corplot_0.9 <- sites_corplot_0.9[sites_corplot_0.9$Var1!=sites_corplot_0.9$Var2,]
length(sites_corplot_0.9$Var1)

# Save correlation_plot as .svg
svglite("figures/corrplot.svg", width = 6.2, height = 6.2)
correlation_plot <- corrplot(sites_corplot, method="circle", type = 'lower', tl.col = 'black', tl.cex = 0.5,
                             cl.ratio = 0.2, tl.srt = 45)
dev.off()

rm(mdata_filtered_sd_spread)
rm(mdata_filtered_sd_spread_na_omit)
rm(mdata_SD_per_site)
rm(correlation_plot)
rm(sites_corplot_0.9)
rm(sites_corplot)

# Combine sites due to high correlation, by taking their mean:
## chr5_71392849 & chr5_71392860 --> chr5_71392849.
## chr5_71409194 & chr5_71409272 & chr5_71409485 --> chr5_71409194.
## chr5_71409976 & chr5_71409990 & chr5_71410038 & chr5_71410055 & chr5_71410068 --> chr5_71409976.
## Adjust in dataframe:
### first select relevant columns and pivot the dataframe to have 1 column per site, then make extra columns for the averages
mdata_filtered_mean_correlated_sites_percentage <- mdata_filtered %>%
  select(-CHROM, -start, -stop, -Me, -Un, -cov, -PCR_prod, -paraphase_SNP, -stdev_percentage_per_site, -percentage_SD_lower_5) %>%
  spread(key = POS_ID, value = percentage) %>%
  mutate(chr5_71392849. = (chr5_71392849 + chr5_71392860) / 2,
         chr5_71409194. = (chr5_71409194 + chr5_71409272 + chr5_71409485) / 3,
         chr5_71409976. = (chr5_71409976 + chr5_71409990 + chr5_71410038 + chr5_71410055 + chr5_71410068) / 5)

### Remove the original columns and pivot the dataframe back
mdata_filtered_mean_correlated_sites_percentage <- mdata_filtered_mean_correlated_sites_percentage %>%
  select(Anonymized_ID, lib_size, SMN_CN_total, lib_size_per_SMN_copy, GQN, SMN1_CN, SMN2_CN, NAIP_CN, mut_859_GC, mut_859_GC_factor, SMN2_RNA, sex, sex_factor, SMA_type, SMA_type_factor, age_at_onset_years, age_at_sampling_years, SMA_type_numeric, CN_minus_type, Concordance, age_group, age_group_factor, chr5_71392849., chr5_71409194., chr5_71409976.) %>%
  gather(chr5_71392849., chr5_71409194., chr5_71409976., key = "POS_ID", value = "percentage") %>%
  left_join(sample_info)

### join with the original table and remove original columns from original table
mdata_filtered_condensed_percentage <- mdata_filtered %>%
  filter(POS_ID != "chr5_71392849" & POS_ID != "chr5_71392860" & POS_ID != "chr5_71409194" & POS_ID != "chr5_71409272" & POS_ID != "chr5_71409485" & POS_ID != "chr5_71409976" & POS_ID != "chr5_71409990" & POS_ID != "chr5_71410038" & POS_ID != "chr5_71410055" & POS_ID != "chr5_71410068") %>%
  full_join(mdata_filtered_mean_correlated_sites_percentage) 
mdata_filtered_condensed_percentage_no_SMN1_859 <- mdata_filtered_condensed_percentage %>%
  filter(SMN1_CN == 0) %>%
  filter(mut_859_GC == "no")
rm(mdata_filtered_mean_correlated_sites_percentage)

# Check if model works at one site for one independent variable
check_model <- run_diff_meth(
  site = "chr5_71393063",
  dependent_variable = "percentage",
  independent_variable = "age_at_sampling_years",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = TRUE,
  method = "lm"
)

summary(check_model$model)
rm(check_model)

# Run model: test if sites are significantly associated with age
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(age_at_sampling_years))$Anonymized_ID)) #number of samples
results_lm_age <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "age_at_sampling_years",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_age <- results_lm_age %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Age at sampling") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-0.6,0.6)) +
  scale_y_continuous(limits = c(0,25)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none")

volcano_age

#number of significant sites
length(unique(filter(results_lm_age, Padj_fdr < 0.01)$site)) 

#significant sites
age_significant_sites <- distinct(select(filter(results_lm_age, Padj_fdr < 0.01), site))

#non-significant sites
age_non_significant_sites <- distinct(select(filter(results_lm_age, Padj_fdr >= 0.01), site)) 

# Save volcano plot for age as .svg 
ggsave(
  "figures/volcano_age.svg",
  plot = volcano_age,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Plot sites that are significantly associated with age
age_significant_sites_plot <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  mutate(significant = case_when(POS_ID %in% filter(results_lm_age, Padj_fdr < 0.01)$site ~ "Yes",
                                 !(POS_ID %in% filter(results_lm_age, Padj_fdr < 0.01)$site) ~ "No")) %>%
  filter(!is.na(age_at_sampling_years)) %>%
  ggplot(aes(x = age_at_sampling_years, y = percentage, color = significant)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5) +
  xlab("Age at sampling (years)") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(values = c("black", "red"), name = "Significant\nassociation\nwith age") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(POS_ID), ncol = 7)
age_significant_sites_plot

# Save sites plot as .svg
ggsave(
  "figures/age_significant_sites.svg",
  plot = age_significant_sites_plot,
  scale = 1,
  width = 16,
  height = 14,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Run model: test if sites are significantly associated with sex (factor)
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(sex))$Anonymized_ID)) #number of samples

results_lm_sex <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "sex_factor",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot for factor sex
volcano_sex <- results_lm_sex %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Sex") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-6,6)) +
  scale_y_continuous(limits = c(0,25)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_sex

# Save volcano_sex as .svg
ggsave(
  "figures/volcano_sex.svg",
  plot = volcano_sex,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with SMN2 copy number
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(SMN2_CN))$Anonymized_ID)) #number of samples

results_lm_SMN2_CN <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMN2_CN",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMN2_CN <- results_lm_SMN2_CN %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2 copy number") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_SMN2_CN

# Save volcano_SMN2_CN as .svg
ggsave(
  "figures/volcano_SMN2_CN.svg",
  plot = volcano_SMN2_CN,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Run model: test if sites are significantly associated with SMA type
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(SMA_type_factor))$Anonymized_ID)) #number of samples
results_lm_SMA_type <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMA_type_factor",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMA_type <- results_lm_SMA_type %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMA type") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_SMA_type

# Save volcano_SMA_type as .svg
ggsave(
  "figures/volcano_SMA_type.svg",
  plot = volcano_SMA_type,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with NAIP copy number
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(NAIP_CN))$Anonymized_ID)) #number of samples
results_lm_NAIP_CN <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "NAIP_CN",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_NAIP_CN <- results_lm_NAIP_CN %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("NAIP copy number") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-3,3)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_NAIP_CN

# Save volcano_NAIP_CN as .svg
ggsave(
  "figures/volcano_NAIP_CN.svg",
  plot = volcano_NAIP_CN,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with SMA type in 3xSMN2 group
mdata_filtered_condensed_percentage_no_SMN1_859_CN3 <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(SMN2_CN == 3)
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859_CN3, !is.na(SMA_type_factor))$Anonymized_ID)) #number of samples
results_lm_SMA_type_CN3 <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859_CN3$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMA_type_factor",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859_CN3,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMA_type_CN3 <- results_lm_SMA_type_CN3 %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMA type (3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-6,6)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_SMA_type_CN3

# Save volcano_SMA_type_CN3 as .svg
ggsave(
  "figures/volcano_SMA_type_CN3.svg",
  plot = volcano_SMA_type_CN3,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Run model: test if sites are significantly associated with SMA type in 4xSMN2 group
mdata_filtered_condensed_percentage_no_SMN1_859_CN4 <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(SMN2_CN == 4)
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859_CN4, !is.na(SMA_type_factor))$Anonymized_ID)) #number of samples
results_lm_SMA_type_CN4 <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859_CN4$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMA_type_factor",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859_CN4,
  returnModel = FALSE,
  method = "lm"
  ) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMA_type_CN4 <- results_lm_SMA_type_CN4 %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMA type (4xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-6,6)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_SMA_type_CN4

# Save volcano_SMA_type_CN4 as .svg
ggsave(
  "figures/volcano_SMA_type_CN4.svg",
  plot = volcano_SMA_type_CN4,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with age at onset
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859, !is.na(age_at_onset_years))$Anonymized_ID)) #number of samples
results_lm_age_at_onset <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "age_at_onset_years",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_age_at_onset <- results_lm_age_at_onset %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Age at onset") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-0.3,0.3)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_age_at_onset

# Save volcano_age_at_onset as .svg
ggsave(
  "figures/volcano_age_at_onset.svg",
  plot = volcano_age_at_onset,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with age at onset in 3xSMN2 group
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859_CN3, !is.na(age_at_onset_years))$Anonymized_ID)) #number of samples

results_lm_age_at_onset_CN3 <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859_CN3$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "age_at_onset_years",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859_CN3,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_age_at_onset_CN3 <- results_lm_age_at_onset_CN3 %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Age at onset (3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_age_at_onset_CN3

# Save volcano_age_at_onset_CN3 as .svg
ggsave(
  "figures/volcano_age_at_onset_CN3.svg",
  plot = volcano_age_at_onset_CN3,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with age at onset in 4xSMN2 group
length(unique(filter(mdata_filtered_condensed_percentage_no_SMN1_859_CN4, !is.na(age_at_onset_years))$Anonymized_ID)) #number of samples
results_lm_age_at_onset_CN4 <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859_CN4$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "age_at_onset_years",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = mdata_filtered_condensed_percentage_no_SMN1_859_CN4,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_age_at_onset_CN4 <- results_lm_age_at_onset_CN4 %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("Age at onset (4xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_age_at_onset_CN4

# Save volcano_age_at_onset_CN4 as .svg
ggsave(
  "figures/volcano_age_at_onset_CN4.svg",
  plot = volcano_age_at_onset_CN4,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Plot SMN2 RNA expression with age
SMN2_RNA_age_graph <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(!is.na(SMN2_RNA)) %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4) %>% #only patients with 3 or 4 copies SMN2, because there is barely any data available ofr the two and five copy group (not enough to do good statistics)
  select(Anonymized_ID, SMN2_CN, age_at_sampling_years, SMN2_RNA, age_group) %>%
  distinct() %>%
  ggplot(aes(x = age_at_sampling_years, y = SMN2_RNA, colour = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_point(size = 1, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  xlab("Age at sampling (years)") +
  ylab("SMN2-FL expression (a.u.)") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), name = "SMN2 copy number") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(SMN2_CN), ncol = 2)

SMN2_RNA_age_graph

# Save SMN2_RNA_age_graph as .svg
ggsave(
  "figures/SMN2_RNA_age.svg",
  plot = SMN2_RNA_age_graph,
  scale = 1,
  width = 12,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Check statistically if there is association of SMN2-FL RNA with age per copy group
## 3xSMN2 copies
SMN2_RNA_age_CN3 <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(!is.na(SMN2_RNA)) %>%
  select(Anonymized_ID, sex, SMN2_CN, age_at_sampling_years, SMN2_RNA) %>%
  filter(SMN2_CN == 3) %>%
  distinct()

model_age_RNA_CN3 <- lm(SMN2_RNA ~ age_at_sampling_years, data = SMN2_RNA_age_CN3)
summary(model_age_RNA_CN3)

## 4xSMN2 copies
SMN2_RNA_age_CN4 <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(!is.na(SMN2_RNA)) %>%
  select(Anonymized_ID, SMN2_CN, age_at_sampling_years, SMN2_RNA) %>%
  filter(SMN2_CN == 4) %>%
  distinct()
model_age_RNA_CN4 <- lm(SMN2_RNA ~ age_at_sampling_years, data = SMN2_RNA_age_CN4)
summary(model_age_RNA_CN4)

# SMN2_FL expression plotted against DNA methylation in sites sifnicicant with age
## Plot sites that are significantly associated with age
mdata_filtered_condensed_percentage_no_SMN1_859_CN3 <- filter(mdata_filtered_condensed_percentage_no_SMN1_859, SMN2_CN == 3)

age_significant_sites_SMN2_FL_RNA <- mdata_filtered_condensed_percentage_no_SMN1_859_CN3 %>%
  mutate(significant = case_when(POS_ID %in% filter(results_lm_age, Padj_fdr < 0.01)$site ~ "Yes",
                                 !(POS_ID %in% filter(results_lm_age, Padj_fdr < 0.01)$site) ~ "No")) %>%
  filter(!is.na(SMN2_RNA)) %>%
  filter(!is.na(age_at_sampling_years)) %>%
  #  filter(POS_ID %in% filter(results_lm_age, Padj_fdr < 0.01)$site) %>%
  ggplot(aes(x = SMN2_RNA, y = percentage, color = significant)) +
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5) +
  xlab("SMN2-FL RNA expression") +
  ylab("Methylation percentage (%)") +
  scale_color_manual(values = c("black", "red"), name = "Significant\nassociation\nwith age") +
  theme_classic(base_size = 7.5) + 
  facet_wrap(vars(POS_ID), ncol = 7)

age_significant_sites_SMN2_FL_RNA

# Save age_significant_sites_SMN2_FL_RNA as svg
ggsave(
  "figures/age_significant_sites_SMN2_FL_RNA.svg",
  plot = age_significant_sites_SMN2_FL_RNA,
  scale = 1,
  width = 16,
  height = 14,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: test if sites are significantly associated with SMN2_FL_RNA (continuous)
results_lm_SMN2_FL_RNA <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_no_SMN1_859_CN3$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMN2_RNA",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"), #
  data = mdata_filtered_condensed_percentage_no_SMN1_859_CN3,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMN2_FL_RNA <- results_lm_SMN2_FL_RNA %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN2-FL RNA expression\n(3xSMN2)") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-15,15)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_SMN2_FL_RNA

# Save volcano_SMN2_FL_RNA as .svg
ggsave(
  "figures/volcano_SMN2_FL_RNA.svg",
  plot = volcano_SMN2_FL_RNA,
  scale = 1,
  width = 4,
  height = 4,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Run model: Test if there is association with 859G>C variant specifically in group with 2 SMN2 copies
results_lm_mut_859 <- purrr::map_df(
  unique(filter(mdata_filtered_condensed_percentage, SMN2_CN == 2)$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "mut_859_GC_factor",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"), 
  data = filter(mdata_filtered_condensed_percentage, SMN2_CN == 2),
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_mut_859 <- results_lm_mut_859 %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("c.859G>C variant") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-60,60)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_mut_859

# Save volcano plot as .svg
ggsave(
  "figures/volcano_mut_859.svg",
  plot = volcano_mut_859,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


# Run model: Test if there is association with SMN1 copy number specifically in group with 2 or 3 SMN2 copies
results_lm_SMN1_CN <- purrr::map_df(
  unique(filter(mdata_filtered_condensed_percentage, SMN2_CN == 2 | SMN2_CN == 3)$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "SMN1_CN",
  ID = "POS_ID",
  covar = c("GQN", "lib_size_per_SMN_copy", "age_at_sampling_years", "sex"),
  data = filter(mdata_filtered_condensed_percentage, SMN2_CN == 2 | SMN2_CN == 3),
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_SMN1_CN <- results_lm_SMN1_CN %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("SMN1 copy number") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-15,15)) + # adjust
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

volcano_SMN1_CN

# Save as .svg
ggsave(
  "figures/volcano_SMN1_CN.svg",
  plot = volcano_SMN1_CN,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

### Treatment response
# Keep only 3 and 4 copy groups because 2 and 5 copy groups are very small
mdata_filtered_condensed_percentage_hfmse <- mdata_filtered_condensed_percentage_no_SMN1_859 %>%
  filter(SMN2_CN == 3 | SMN2_CN == 4)

# Median age in dHFMSE group
mdata_filtered_condensed_percentage_hfmse %>%
  filter(!is.na(dHFMSE)) %>%
  select(Anonymized_ID, age_at_sampling_years) %>%
  distinct() %>%
  filter(!is.na(age_at_sampling_years)) %>%
  summarise(median_age = median(age_at_sampling_years),
            min_age = min(age_at_sampling_years),
            max_age = max(age_at_sampling_years))

# Determine and plot baseline scores
HFMSE_baseline_scores <- mdata_filtered_condensed_percentage_hfmse %>%
  select(Anonymized_ID, HFMSE_baseline, dHFMSE, SMN2_CN, HFMSE_response_group) %>%
  filter(!is.na(dHFMSE)) %>%
  distinct() %>%
  ggplot(aes(x = HFMSE_baseline, y = dHFMSE, colour = as.factor(SMN2_CN), fill = as.factor(SMN2_CN))) +
  geom_hline(yintercept = -2.5, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  xlab("Baseline HFMSE score") +
  ylab("dHFMSE") +
  labs(fill = "SMN2 CN", color = "SMN2 CN") +
  scale_colour_manual(values = colour_palette) +
  scale_fill_manual(values = colour_palette) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )

HFMSE_baseline_scores

# Save HFMSE baseline scores plot as svg
ggsave(
  "figures/HFMSE_baseline_scores.svg",
  plot = HFMSE_baseline_scores,
  scale = 1,
  width = 7.5,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)

# Run model: test if sites are significantly associated with HFMSE response
length(unique(filter(mdata_filtered_condensed_percentage_hfmse, !is.na(dHFMSE))$Anonymized_ID))

results_lm_HFMSE <- purrr::map_df(
  unique(mdata_filtered_condensed_percentage_hfmse$POS_ID),
  .f = run_diff_meth,
  dependent_variable = "percentage",
  independent_variable = "dHFMSE",
  ID = "POS_ID",
  covar = c("GQN", "age_at_sampling_years", "lib_size_per_SMN_copy", "sex"), 
  data = mdata_filtered_condensed_percentage_hfmse,
  returnModel = FALSE,
  method = "lm"
) %>%
  mutate(Padj_fdr = p.adjust(P,"fdr"))

# Make volcano plot
volcano_HFMSE <- results_lm_HFMSE %>%
  mutate(Padj_fdr_sig = case_when(Padj_fdr < 0.01 ~ "yes",
                                  Padj_fdr >= 0.01 ~ "no")) %>%
  ggplot(aes(y = -log10(Padj_fdr), x = b, color = Padj_fdr_sig)) + 
  geom_point(size = 1, stroke = 0, alpha = 0.5) +
  geom_hline(yintercept =2.00000, color = "black", linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  ggtitle("dHFMSE") +
  xlab("Estimate") + 
  ylab(expression(paste(-log[10], "(P"[adj],")"))) +
  scale_x_continuous(limits = c(-1.1,1.1)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic(base_size = 7.5) + 
  theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0.2)) +
  theme(legend.position = "none"
  )
volcano_HFMSE

# Save HFMSE volcano plot as .svg
ggsave(
  "figures/volcano_HFMSE.svg",
  plot = volcano_HFMSE,
  scale = 1,
  width = 4,
  height = 3.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#write results in table
write_tsv(results_lm_sex, "results_tables/3C_results_lm_sex.txt")
write_tsv(results_lm_age, "results_tables/3D_results_lm_age.txt")
write_tsv(results_lm_SMN2_CN, "results_tables/4C_results_lm_SMN2_CN.txt")
write_tsv(results_lm_SMA_type, "results_tables/4D_results_lm_SMA_type.txt")
write_tsv(results_lm_HFMSE, "results_tables/4G_results_lm_HFMSE.txt")
write_tsv(results_lm_SMN2_FL_RNA, "results_tables/S7C_results_lm_SMN2_FL_RNA.txt")
write_tsv(results_lm_SMN1_CN, "results_tables/S8A_results_lm_SMN1_CN.txt")
write_tsv(results_lm_mut_859, "results_tables/S8B_results_lm_mut_859.txt")
write_tsv(results_lm_NAIP_CN, "results_tables/S9A_results_lm_NAIP_CN.txt")
write_tsv(results_lm_age_at_onset, "results_tables/S9B_results_lm_age_at_onset.txt")
write_tsv(results_lm_SMA_type_CN3, "results_tables/S9C_results_lm_SMA_type_CN3.txt")
write_tsv(results_lm_SMA_type_CN4, "results_tables/S9D_results_lm_SMA_type_CN4.txt")
write_tsv(results_lm_age_at_onset_CN3, "results_tables/S9E_results_lm_age_at_onset_CN3.txt")
write_tsv(results_lm_age_at_onset_CN4, "results_tables/S9F_results_lm_age_at_onset_CN4.txt")

#clean up variables
rm(results_lm_sex)
rm(results_lm_age)
rm(results_lm_SMN2_CN)
rm(results_lm_SMA_type)
rm(results_lm_HFMSE)
rm(results_lm_SMN2_FL_RNA)
rm(results_lm_NAIP_CN)
rm(results_lm_age_at_onset)
rm(results_lm_SMA_type_CN3)
rm(results_lm_SMA_type_CN4)
rm(results_lm_age_at_onset_CN3)
rm(results_lm_age_at_onset_CN4)
rm(results_lm_SMN1_CN)
rm(results_lm_mut_859)

rm(volcano_sex)
rm(volcano_age)
rm(volcano_SMN2_CN)
rm(volcano_SMA_type)
rm(volcano_HFMSE)
rm(volcano_SMN2_FL_RNA)
rm(volcano_NAIP_CN)
rm(volcano_age_at_onset)
rm(volcano_SMA_type_CN3)
rm(volcano_SMA_type_CN4)
rm(volcano_age_at_onset_CN3)
rm(volcano_age_at_onset_CN4)
rm(volcano_SMN1_CN)
rm(volcano_mut_859)

rm(age_significant_sites)
rm(age_non_significant_sites)
rm(age_significant_sites_plot)
rm(age_significant_sites_SMN2_FL_RNA)
rm(HFMSE_baseline_scores)
rm(SMN2_RNA_age_CN3)
rm(SMN2_RNA_age_CN4)
rm(SMN2_RNA_age_graph)
rm(model_age_RNA_CN3)
rm(model_age_RNA_CN4)

