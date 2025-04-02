# script for filtering microbiome data

library(tidyverse)
library(dplyr)
library(magrittr)
library(Biostrings)
library(DECIPHER)

rm(list=ls())


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


# combining asv data and SM data
data_asv_f <- data_asv %>% 
  dplyr::select(host_ID, unfiltered_reads, contains("ASV"))

dat <- left_join(data_sm, data_asv_f, by="host_ID") 


########################################
# start the filtering

##### filter 1
# filtering non-rattus host species
dat1 <- dat %>% 
  filter(host_species == "Rattus rattus") %>% 
  select_if(~ any(. != 0))  # removing all ASVs not belonging to rattus (0 in all samples)
# here I filtered 12357 ASVs
# we start with 10352 ASVs

##### filter 2
# removing ASVs based on taxonomy

# adding taxonomy
tax <- read_delim("data/data_raw/data_microbiome/ASVs_taxonomy_new.tsv") %>% 
  dplyr::rename(asv_ID = ASV)

# setting the not allowed taxonomy:
# not bacteria, chloroplast or mitochondria
tax_exclude <- tax %>% 
  filter(asv_ID %in% colnames(dat1)) %>% 
  filter(Kingdom != "Bacteria" | Order == "Chloroplast" | Family == "Mitochondria" | is.na(Kingdom))
# here I filtered 68 ASVs

dat2 <- dat1 %>% 
  select(-all_of(tax_exclude$asv_ID))


##### filter 3
# removing ASVs with very low relative read abundance in each sample

asv_rel_reads_th <- 0.001 # threshold of 0.1%
dat3 <- dat2 %>%
  mutate(across(starts_with("ASV"),~ ./unfiltered_reads)) %>% 
  mutate(across(starts_with("ASV"), ~ifelse(.<asv_rel_reads_th,0,.))) %>% 
  mutate(across(starts_with("ASV"),~ .*unfiltered_reads)) %>%
  select_if(~ any(. != 0))
# here I filtered 3255 ASVs



##### filter 4
# removing low occurrence ASVs

# transforming to long format
dat3_long <- dat3 %>% 
  gather("asv_ID", "reads", starts_with("ASV")) %>% 
  filter(reads>0)

# number of hosts
n_host <- length(unique(dat3_long$host_ID)) 

# ASVs prevalence across all villages and grids
asv_occur <- dat3_long %>% 
  count(asv_ID) %>% 
  mutate(n_host = n_host) %>% 
  mutate(host_p = n/n_host)

# plotting
asv_occur %>% 
  ggplot(aes(x=n)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
  labs(x="No. of occurrences", y="Count")

# finding the best filter
asv_unique <- NULL
for (i in seq(0.01,0.8, by=0.005)) {
  asv_unique_i <- asv_occur %>% 
    filter(host_p>i) %>% 
    summarise(asv_n = n_distinct(asv_ID)) %>% 
    mutate(i=i)
  
  asv_unique <- rbind(asv_unique, asv_unique_i)
}
#pdf(file = 'results/figure_S1.pdf')
asv_unique %>% 
  ggplot(aes(x=i, y=asv_n/1951)) + 
  geom_line() +
  geom_point() +
  geom_vline(xintercept=0.2, linetype='dashed', color="gray") +
  geom_hline(yintercept=103/1951, linetype='dashed', color="gray") +
  #scale_y_continuous(limits = c(0, 0.35)) +
  #scale_x_continuous(limits = c(0, 20)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 15), panel.grid = element_blank(), panel.border = element_rect(color = "black")) +
  labs(x="ASV Prevalence", y="Prop. of ASVs")
#dev.off()

# filtering the ASVs
# threshold of 1% prevalence
asv_occur_th <- 0.01

dat4 <- asv_occur %>% 
  filter(host_p > asv_occur_th) %>% 
  select(asv_ID, host_p) %>% 
  left_join(dat3_long, by=c("asv_ID"))
# here I filtered 5078 ASVs



###
# calculating the new final total reads per sample
host_total_reads <- dat4 %>% 
  group_by(host_ID) %>% 
  summarise(total_reads = sum(reads))

dat4 %<>% left_join(host_total_reads, by="host_ID") %>% 
  select(-unfiltered_reads) %>% 
  mutate(reads_p = reads/total_reads)


# how many ASVs?
length(unique(dat4$asv_ID)) # 1951 ASVs
dat4 %>% group_by(village) %>% summarise(n_distinct(asv_ID))
# host richness
host_richness <- dat4 %>% group_by(host_ID) %>% summarise(n=n_distinct(asv_ID))
hist(host_richness$n)


##### filter 5
# removing samples (hosts) with low total reads

# plotting total reads distribution
dat4 %>% 
  distinct(host_ID, total_reads) %>% 
  ggplot(aes(x=total_reads)) +
  geom_histogram(binwidth = 1000) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
  labs(x="Total Reads", y="Count")

# removing samples with less than 5000 total reads
total_reads_th <- 5000
dat5 <- dat4 %>% 
  filter(total_reads > total_reads_th) %>% 
  mutate(reads = reads_p) %>% 
  select(-reads_p)

# how many hosts we lose?
length(unique(dat4$host_ID)) - length(unique(dat5$host_ID)) # 17 hosts
# how many ASVs?
length(unique(dat5$asv_ID)) # 1951 ASVs


# plotting prevalence distribution
n_host_final <- length(unique(dat5$host_ID))

dat5 %>% 
  count(asv_ID) %>% 
  mutate(asv_p = n/n_host_final) %>% 
  ggplot(aes(x=asv_p, y=..ncount..)) + 
  geom_histogram() +
  geom_vline(xintercept=c(0.2), linetype='dashed', color="gray") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 15), panel.grid = element_blank(), panel.border = element_rect(color = "black")) +
  labs(x="ASV Prevalence", y="No. of ASVs")


# saving the data
write_csv(dat5, "data/data_processed/data_asv_rra0.001_p0.01_th5000_all.csv")


