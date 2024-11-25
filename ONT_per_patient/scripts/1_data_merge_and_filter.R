# Load required libraries for "1_data_merge_and_filter.R"
library(tidyverse)
library(stringr)

# Set working directory
## Set working directory to folder ONT_per_patient, to execute all commands without adjusting the paths of the files

# Load in sample info
sample_info <- read_delim("data/sample_info_ONT_complete_anonymized.txt", delim = "\t")

# Load in merged modbam2bed files with extra column for sample
bed_full <- read_delim("data/modbam2bed_ONT_per_patient_anonymized.bed", delim = "\t", col_names = FALSE)

# Adjust column names of the bed file
colnames(bed_full) <- c("CHROM", "start", "POS", "modbase", "score", "strand",
                        "remove1", "remove2", "remove3",
                        "coverage", "modpercentage", "Ncanonical", "Nmodified",
                        "Nfilter", "sample")

bed_full <- bed_full %>%
  select(CHROM, start, POS, strand, coverage, modpercentage, Ncanonical, Nmodified, sample) %>%
  mutate(POS = case_when(strand == "-" ~ POS - 1,
                         strand == "+" ~ POS)) %>%
  #change coordinate of the - strand rows, so that they correspond to the + strand and can be merged later
  mutate(start = POS)

# Merging the "-" and "+" strand
bed_strands_merged <- bed_full %>%
  mutate(CHROM = as.numeric(str_remove(CHROM, "chr"))) %>%
  group_by(sample, CHROM, start, POS) %>%
  summarise(Nmod = sum(Nmodified),
            Ncan = sum(Ncanonical),
            cov_tot = sum(coverage)) %>%
  mutate(cov_called = Ncan + Nmod,
         percentage = Nmod / cov_called * 100) %>%
  # attention: this percentage is calculated over all bases that have a methylation call.
  # If you want to include the filtered bases in the total, divide by cov_tot
  mutate(Anonymized_ID = str_extract(sample, "SMA_\\d{2,3}"),
         tissue = case_when(grepl("fib", sample) ~ "fib",
                            grepl("blood", sample) ~ "blood")) %>%
  mutate(Anonymized_ID_tissue = paste(Anonymized_ID, tissue, sep = "_")) %>%
  ungroup() %>%
  left_join(sample_info, by = "Anonymized_ID") %>%
  mutate(age_at_sampling = case_when(tissue == "fib" ~ age_at_biopsy_years,
                                     tissue == "blood" ~ age_at_EDTA_DNA_years)) %>%
  mutate(age_group = case_when(age_at_sampling < 18 ~ "pediatric",
                               age_at_sampling >= 18 ~ "adult"))

# Calculate library size per sample with coverage (total)
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
