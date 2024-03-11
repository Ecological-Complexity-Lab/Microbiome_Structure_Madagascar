# script for filtering microbiome data

library(tidyverse)
library(dplyr)

# reading small mammals data
data_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")

# reading ASVs raw data
data_asv <- read_csv("data/data_raw/data_microbiome/merged_full_sample_table.csv") 
data_asv %<>% filter(sample_type == "SAMPLE") %>% 
  filter(!(grepl("D",Sample_Name)))

# extracting samples IDs (taking only numbers due to unmatching formats)
data_mammals %<>% mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id)))
data_asv %<>% mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name)))

# matching small-mammals (SM) IDs in the two data files and taking only SM with microbes data
data_sm <- semi_join(data_mammals, data_asv, by="host_ID") %>% 
  dplyr::select(host_ID, final_identification, village, habitat_type, season) %>%
  dplyr::rename(host_species = final_identification, grid = habitat_type) %>% 
  mutate(season = factor(season, levels = c("1","2","3"))) %>% 
  mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth","agriculture","flooded_rice","agroforest","village")))

# combining asv data and SM data
data_asv_f <- data_asv %>% 
  dplyr::select(host_ID, unfiltered_reads, contains("ASV"))

dat <- left_join(data_sm, data_asv_f, by="host_ID") 

### filtering out ASVs with less than threshold of relative abundance
rel_reads_threshold <- 0.01
dat_filtered_rel_abundance <- dat %>%
  mutate(across(starts_with("ASV"),~ ./unfiltered_reads)) %>% 
  mutate(across(starts_with("ASV"), ~ifelse(.<rel_reads_threshold,0,.)))

# calculate updated total reads
dat_filtered_rel_abundance %<>% mutate(across(starts_with("ASV"), ~.*unfiltered_reads)) %>% 
  mutate(total_reads = select(.,starts_with("ASV")) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(across(starts_with("ASV"), ~./total_reads))

# filtering out samples with less than threshold of reads
reads_threshold <- 1000
dat_filtered_threshold <- dat_filtered_rel_abundance %>% 
  filter(total_reads >= reads_threshold)

# transforming to long format
data_asv_long_format <- dat_filtered_threshold %>% 
  relocate(total_reads, .before = unfiltered_reads) %>% 
  select(-unfiltered_reads) %>% 
  gather("asv_ID", "reads", starts_with("ASV")) %>% 
  filter(reads>0)

write_csv(data_asv_long_format, "data_processed/three_villages/data_asv_rra0.01_th1000.csv")
