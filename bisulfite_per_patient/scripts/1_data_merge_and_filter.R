# Load the required libraries
library(tidyverse)

# Set working directory to the bisulfite_analysis_R_20241119 folder to use all the paths without changing them

# Load the sample info table
sample_info <- read_delim("data/sample_info_complete.txt", delim = "\t")

# Assign SMA_ID as rownames and check column types
rownames(sample_info) <- sample_info$SMA_ID
lapply(sample_info, class)

# Make a list of the .cov.gz files that contain the methylation data
file_list <- list.files(path = "./data/methylation_cov_gz_files",
                        pattern = "*.cov.gz", full.names = TRUE)

# Extract sample names from file_list
sample_names <- gsub("\\_1_val_1_bismark_bt2_pe.bismark.cov.gz$", "", basename(file_list))
sample_names

# Make list containing all data; 
data_list <- lapply(seq_along(file_list), function(i) {
  # Add column with sample_name to dataframe
  data <- read.delim(file_list[i], sep = "\t", header = FALSE)
  data$sample <- sample_names[i]
  return(data)
})

# Combine all dataframes with methylation data into a single dataframe
combined_data <- bind_rows(data_list)

## Rename column names
colnames(combined_data) <- c("CHROM", "start", "stop", "percentage", "Me", "Un", "SMA_ID")

## Compute library size for each sample
lib_size <- combined_data %>%
  group_by(SMA_ID) %>%
  summarise(lib_size = sum(Me+Un))

# Add library size to the dataframe with the methylation data
combined_data <- combined_data %>%
  left_join(lib_size, by = "SMA_ID")

# Add columns to the dataframe and merge with the sample_info table
combined_data <- combined_data %>%
  mutate(POS_ID = paste0(CHROM, "_", stop, sep = "")) %>% #make a character type column for genomic position
  mutate(cov = Me + Un) %>% #calculate total coverage
  mutate(PCR_prod = case_when(between(stop, 71380875, 71381350) ~ "SMN_prom_0.2",
                              between(stop, 71381211, 71381698) ~ "SMN_prom_1.3",
                              between(stop, 71381682, 71382233) ~ "SMN_prom_2.2",
                              between(stop, 71382327, 71382824) ~ "SMN_prom_3.5", #has 6 sites that overlap with SMN_prom_4.2, they are included in this amplicon in the annotation
                              between(stop, 71382825, 71383133) ~ "SMN_prom_4.2", #has 6 sites that overlap with SMN_prom_3.5, they are not included in this amplicon in the annotation
                              between(stop, 71385508, 71386053) ~ "SMN_intron1_1.1",
                              between(stop, 71387871, 71388345) ~ "SMN_lncRNA_1.3",
                              between(stop, 71388354, 71388807) ~ "SMN_lncRNA_2.1", #has 4 sites that overlap with SMN_lncRNA_3a.1, they are included in this amplicon in the annotation
                              between(stop, 71388808, 71388975) ~ "SMN_lncRNA_3a.1", #has 4 sites that overlap with SMN_lncRNA_2.1, they are not included in this amplicon in the annotation
                              between(stop, 71388956, 71389455) ~ "SMN_lncRNA_3b.3",
                              between(stop, 71389434, 71389903) ~ "SMN_lncRNA_4.3",
                              between(stop, 71390145, 71390644) ~ "SMN_lncRNA_5.2",
                              between(stop, 71392682, 71393138) ~ "SMN_intron1_2.1",
                              between(stop, 71394810, 71395294) ~ "SMN_intron1_3.1",
                              between(stop, 71397422, 71397921) ~ "SMN_intron2a_1.1",
                              between(stop, 71406895, 71407431) ~ "SMN_intron6_2.1",
                              between(stop, 71409095, 71409544) ~ "SMN_3-UTR_1.2",
                              between(stop, 71409656, 71410119) ~ "SMN_3-UTR_2.1")) %>%
  mutate(paraphase_SNP = case_when(stop %in% c(71380908, 71381869, 71381906, 71382006,
                                               71382644, 71385830, 71385998, 71388149,
                                               71388417, 71390661, 71390689, 71394993,
                                               71397505, 71406980, 71407128, 71407136,
                                               71407280, 71409298, 71410048) ~ "yes")) %>% #CpG site locations that are a known SNV site
  left_join(sample_info) %>%
  mutate(lib_size_per_SMN_copy = lib_size / SMN_CN_total) #add library size per SMN copy, because the more SMN copies a sample has, the higher the read depth likely is.

# Filter combined data as described in figure 
combined_data_filtered <- combined_data %>%
  filter(!is.na(PCR_prod)) %>% #filter for only sites on SMN PCR product sites --> 385 sites left
  filter(cov >= 100) %>% #filter to only contain CpG sites on the SMN gene (with flanks) and with >= 100x coverage
  filter(is.na(paraphase_SNP))  #filter out paraphase SNPs

# check metrics
## total number of sites: 55142
length(unique(combined_data$POS_ID)) 

## number of on-target CpG sites with cov >=100: 129
length(unique(filter(filter(combined_data, !is.na(PCR_prod)), cov >= 100)$POS_ID)) 

## number of sites after filtering out paraphase SNPs: 121
length(unique(filter(filter(filter(combined_data, !is.na(PCR_prod)), cov >= 100), is.na(paraphase_SNP))$POS_ID)) 

# clean up variables
rm(file_list)
rm(sample_names)
rm(data_list)
rm(lib_size)
