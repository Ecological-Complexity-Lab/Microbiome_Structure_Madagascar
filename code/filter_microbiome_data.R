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
  dplyr::select(host_ID, field_identification, village, habitat_type, season) %>%
  dplyr::rename(host_species = field_identification, grid = habitat_type) %>% 
  mutate(season = factor(season, levels = c("1","2","3"))) %>% 
  mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth","agriculture","flooded_rice","agroforest","village")))

# saving attributes only for relevant small mammals
attr_SM <- data_sm %>% 
  left_join(data_mammals %>% select(host_ID,latitude,longitude,elevation.obs,sex,mass,age_repro,age_dental), by="host_ID")
# writing to a csv file
write_csv(attr_SM, "data/data_processed/small_mammals/small_mammals_attributes.csv")


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

#write_csv(data_asv_long_format, "data/data_processed/microbiome/data_asv_rra0.01_th1000.csv")


#####
# adding taxonomy
tax <- read_delim("../data/data_raw/data_microbiome/ASVs_taxonomy_new.tsv") %>% 
  dplyr::rename(asv_ID = ASV)

# setting the not allowed taxonomy:
# not bacteria, chloroplast or mitochondria
tax_exclude <- tax %>% 
  filter(Kingdom != "Bacteria" | Order == "Chloroplast" | Family == "Mitochondria" | is.na(Kingdom))

# finding phylogenetic distance
# reading the phylogenetic tree
best_tree <- readRDS(file = "results/phylo_tree_0.01_2.rds")
phylo_tree <- best_tree$tree
# ASVs phylogenetic distance
asv_distance <- ape::cophenetic.phylo(phylo_tree)

# mean distance for each ASV
mean_phylo_dist <- rowMeans(asv_distance)
hist(mean_phylo_dist)
a <- names(mean_phylo_dist[mean_phylo_dist>3])
b <- tax %>% filter(!(asv_ID %in% tax_exclude$asv_ID) & asv_ID %in% a)
# Not many ASVs seem very different from all the rest (dist>4). Only ASV_4819 is not in one of the categories of the not allowed taxonomy. 
# So I remove it from the analysis. 

data_asv_long_format <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv")

data_asv_long_format_clean <- data_asv_long_format %>% 
  filter(!(asv_ID %in% tax_exclude$asv_ID)) %>% 
  filter(asv_ID != "ASV_4819")


write_csv(data_asv_long_format_clean, "data/data_processed/microbiome/data_asv_rra0.01_th1000_clean.csv")
