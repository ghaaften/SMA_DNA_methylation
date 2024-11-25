# Load required libraries for "1_process_bed_file.R"
library(tidyverse)
library(stringr)

# Set working directory to ONT_per_haplotype folder to run all commands without changing file paths

# Load in sample info
sample_info <- read_delim("data/sample_info_ONT_complete_anonymized.txt", delim = "\t")

# Load in SMN downstream environment info
hybrid_type_downstream_env <- read_delim("data/downstream_env_per_haplotype_anonymized.txt", delim = "\t")

# Load in merged modbam2bed files with extra column for sample
bed_full <- read_delim("data/modbam2bed_ONT_per_haplotype_anonymized.bed", delim = "\t", col_names = FALSE)

# Adjust column names of bed_full df
colnames(bed_full) <- c("CHROM", "start", "POS", "modbase", "score", "strand",
                        "remove1", "remove2", "remove3",
                        "coverage", "modpercentage", "Ncanonical", "Nmodified",
                        "Nfilter", "sample_hap")

bed_full <- bed_full %>%
  select(CHROM, start, POS, strand, coverage, modpercentage, Ncanonical, Nmodified, sample_hap) %>%
  mutate(POS = case_when(strand == "-" ~ POS - 1,
                         strand == "+" ~ POS)) %>% 
  #change coordinate of the - strand rows, so that they correspond to the + strand and can be merged later
  mutate(start = POS)

# # Merging the "-" and "+" strand, add clinical info and including all relevant columns
bed_strands_merged <- bed_full %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "chr"))) %>%
  group_by(sample_hap, CHROM, start, POS) %>%
  summarise(Nmod = sum(Nmodified),
            Ncan = sum(Ncanonical),
            cov_tot = sum(coverage)) %>%
  mutate(cov_called = Ncan + Nmod,
         percentage = Nmod / cov_called * 100) %>% 
         # attention: this percentage is calculated over all bases that have a methylation call. 
         # If you want to include the filtered bases in the total, divide by cov_tot
  mutate(percentage = case_when(cov_called <4 ~ NA,
                                    cov_called >=4 ~ percentage)) %>% 
         # change percentage to NA if cov_called is lower than 4
  mutate(Anonymized_ID = str_extract(sample_hap, "SMA_\\d{2,3}"),
         tissue = case_when(grepl("fib", sample_hap) ~ "fib",
                            grepl("blood", sample_hap) ~ "blood"),
         hap = case_when(grepl("hap1", sample_hap) ~ "hap1",
                         grepl("hap2", sample_hap) ~ "hap2",
                         grepl("hap3", sample_hap) ~ "hap3",
                         grepl("hap4", sample_hap) ~ "hap4",
                         grepl("hap5", sample_hap) ~ "hap5")) %>%
  mutate(Anonymized_ID_tissue_hap = paste(Anonymized_ID, tissue, hap, sep = "_")) %>%
  ungroup() %>%
  left_join(sample_info, by = "Anonymized_ID") %>%
  left_join(hybrid_type_downstream_env, by = c("Anonymized_ID", "tissue", "hap")) %>%
  mutate(age_at_sampling = case_when(tissue == "fib" ~ age_at_biopsy_years,
                                      tissue == "blood" ~ age_at_EDTA_DNA_years)) %>%
  mutate(age_group = case_when(age_at_sampling < 18 ~ "pediatric",
                               age_at_sampling >= 18 ~ "adult"))

# Create lib_size df: calculate library size per sample with cov_tot
lib_size <- bed_strands_merged %>%
  group_by(Anonymized_ID) %>%
  summarise(lib_size = sum(cov_tot))

# Add library size per SMN copy and other annotation per CpG site
bed_strands_merged <- bed_strands_merged %>%
  left_join(lib_size) %>%
  mutate(SMN_CN_total = SMN1_CN + SMN2_CN) %>%
  mutate(lib_size_per_SMN_copy = lib_size / SMN_CN_total) %>%
  mutate(POS_ID = paste("chr5_", POS, sep = "")) %>%
  mutate(PCR_prod = case_when(between(POS, 71380875, 71381350) ~ "SMN_prom_0.2",
                              between(POS, 71381211, 71381698) ~ "SMN_prom_1.3",
                              between(POS, 71381682, 71382233) ~ "SMN_prom_2.2",
                              between(POS, 71382327, 71382824) ~ "SMN_prom_3.5", #has 6 sites that overlap with SMN_prom_4.2, they are included in this amplicon in the annotation
                              between(POS, 71382825, 71383133) ~ "SMN_prom_4.2", #has 6 sites that overlap with SMN_prom_3.5, they are not included in this amplicon in the annotation
                              between(POS, 71385508, 71386053) ~ "SMN_intron1_1.1",
                              between(POS, 71387871, 71388345) ~ "SMN_lncRNA_1.3",
                              between(POS, 71388354, 71388807) ~ "SMN_lncRNA_2.1", #has 4 sites that overlap with SMN_lncRNA_3a.1, they are included in this amplicon in the annotation
                              between(POS, 71388808, 71388975) ~ "SMN_lncRNA_3a.1", #has 4 sites that overlap with SMN_lncRNA_2.1, they are not included in this amplicon in the annotation
                              between(POS, 71388956, 71389455) ~ "SMN_lncRNA_3b.3",
                              between(POS, 71389434, 71389903) ~ "SMN_lncRNA_4.3",
                              between(POS, 71390145, 71390644) ~ "SMN_lncRNA_5.2",
                              between(POS, 71392682, 71393138) ~ "SMN_intron1_2.1",
                              between(POS, 71394810, 71395294) ~ "SMN_intron1_3.1",
                              between(POS, 71397422, 71397921) ~ "SMN_intron2a_1.1",
                              between(POS, 71406895, 71407431) ~ "SMN_intron6_2.1",
                              between(POS, 71409095, 71409544) ~ "SMN_3-UTR_1.2",
                              between(POS, 71409656, 71410119) ~ "SMN_3-UTR_2.1")) %>%
  mutate(Bisulfite_amplicon = case_when(!is.na(PCR_prod) ~ "yes",
                                        is.na(PCR_prod) ~ "no"))

# clean up environment
rm(bed_full)
rm(lib_size)

