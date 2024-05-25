# script for filtering microbiome data

library(tidyverse)
library(dplyr)
library(magrittr)
library(iNEXT)
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

# saving attributes only for relevant small mammals
attr_SM <- data_sm %>% 
  left_join(data_mammals %>% select(host_ID,latitude,longitude,elevation.obs,sex,mass,age_repro,age_dental), by="host_ID")
# writing to a csv file
#write_csv(attr_SM, "data/data_processed/small_mammals/small_mammals_attributes.csv")


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
asv_rel_reads_th <- 0.001
dat3 <- dat2 %>%
  mutate(across(starts_with("ASV"),~ ./unfiltered_reads)) %>% 
  mutate(across(starts_with("ASV"), ~ifelse(.<asv_rel_reads_th,0,.))) %>% 
  mutate(across(starts_with("ASV"),~ .*unfiltered_reads)) %>%
  select_if(~ any(. != 0))
# here I filtered 8546 ASVs


# asv_total_reads_th <- 0.01
# asv_total_reads <- colSums(dat1 %>% select(starts_with("ASV")))
# asv_total_reads_low <- asv_total_reads[asv_total_reads<asv_total_reads_th]
# here (th=100) I filtered 5098 ASVs

# dat2 <- dat1 %>% 
#   select(!names(asv_total_reads_low))



##### filter 4
# removing low occurrence ASVs

# transforming to long format
dat3_long <- dat3 %>% 
  gather("asv_ID", "reads", starts_with("ASV")) %>% 
  filter(reads>0)

# number of hosts in each village
n_host_village <- dat3_long %>% 
  group_by(village) %>% summarise(n_host = n_distinct(host_ID))

# ASVs occurrence by village
asv_occur_village <- dat3_long %>% 
  count(village, asv_ID) %>% 
  left_join(n_host_village, by="village") %>% 
  mutate(host_p = n/n_host)

# plotting
asv_occur_village %>% 
  ggplot(aes(x=n)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~village)+
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
  labs(x="No. of occurrences", y="Count")

# finding the best filter
asv_unique <- NULL
for (i in seq(0,1, by=0.01)) {
  asv_unique_i <- asv_occur_village %>% 
    filter(host_p>i) %>% 
    group_by(village) %>% 
    summarise(asv_n = n_distinct(asv_ID)) %>% 
    mutate(i=i)
  
  asv_unique <- rbind(asv_unique, asv_unique_i)
}

asv_unique %>% 
  ggplot(aes(x=i, y=asv_n)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~village) +
  #scale_y_continuous(limits = c(0, 0.35)) +
  #scale_x_continuous(limits = c(0, 20)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 15), panel.grid = element_blank(), panel.border = element_rect(color = "black")) +
  labs(x="ASV Prevalence", y="No. of ASVs")

# filtering the ASVs
asv_occur_th <- 0.02

dat4 <- asv_occur_village %>% 
  filter(host_p > asv_occur_th) %>% 
  select(village, asv_ID) %>% 
  left_join(dat3_long, by=c("village","asv_ID"))



##### filter 5
# removing very phylogeneticly different ASVs

# finding phylogenetic distance
# reading the full dna sequences (fasta file)
seq_fa <- read.FASTA(file="data/data_raw/data_microbiome/ASV_merged_full.fa")
asv_names <- tibble(asv_ID = unique(dat4$asv_ID))
asv_names %<>% mutate(ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", asv_ID))) %>% arrange(ID)
seq_fa2 <- seq_fa[names(seq_fa) %in% asv_names$asv_ID]
names(seq_fa2) <- asv_names$asv_ID
#write.FASTA(seq_fa2, file = "data/data_raw/data_microbiome/ASV_filtered_new.fa") # writing the filtered fasta file

# aligning the sequences
seq_aligned <- readDNAStringSet("data/data_raw/data_microbiome/ASV_filtered_new.fa")
aligned <- DECIPHER::AlignSeqs(seq_aligned)
seq_aligned2 <- as.DNAbin(aligned)
asv_distance <- as.matrix(ape::dist.dna(seq_aligned2, model = "TN93"))


# mean distance for each ASV
mean_phylo_dist <- rowMeans(asv_distance)
hist(mean_phylo_dist)
quantile(mean_phylo_dist,0.99)
a <- names(mean_phylo_dist[mean_phylo_dist>0.35])
b <- tax %>% filter(!(asv_ID %in% tax_exclude$asv_ID) & asv_ID %in% a)
# Not many ASVs seem very different from all the rest (dist>4). Only ASV_4819 is not in one of the categories of the not allowed taxonomy. 
# So I remove it from the analysis. 

dat3 %<>% select(-ASV_4819)


###
# calculating the new final total reads per sample
host_total_reads <- dat4 %>% 
  group_by(host_ID) %>% 
  summarise(total_reads = sum(reads))

dat4 %<>% left_join(host_total_reads, by="host_ID") %>% 
  select(-unfiltered_reads) %>% 
  mutate(reads_p = reads/total_reads)


# how many ASVs?
length(unique(dat4$asv_ID))
dat4 %>% group_by(village) %>% summarise(n_distinct(asv_ID))
# host richness
host_richness <- dat4 %>% group_by(host_ID) %>% summarise(n=n_distinct(asv_ID))
hist(host_richness$n)


##### filter 6
# removing samples with low total reads

dat4 %>% 
  distinct(host_ID, total_reads) %>% 
ggplot(aes(x=total_reads)) +
  geom_histogram(binwidth = 1000) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
  labs(x="Total Reaads", y="Count")

# accumulation curve
# transforming to matrix
dat4_mat <- dat4 %>% 
  filter(total_reads < 10000) %>% 
  select(host_ID, asv_ID, reads) %>% 
  spread(asv_ID, reads, fill = 0) %>% 
  column_to_rownames("host_ID") %>% 
  select_if(~ any(. != 0)) %>% 
  as.matrix()

 dat4_list <- split(dat4_mat, seq(nrow(dat4_mat)))

#calculating accumulation curve
accum_curve2 <- iNEXT(dat4_list)
accum_curve <- iNEXT(dat4_list) #long
plot(accum_curve2)

# removing samples with less than 5000 total reads
total_reads_th <- 5000
dat5 <- dat4 %>% 
  filter(total_reads > total_reads_th) %>% 
  mutate(reads = reads_p) %>% 
  select(-reads_p)

# how many hosts we lose?
length(unique(dat4$host_ID)) - length(unique(dat5$host_ID))

# saving the data
write_csv(dat5, "data/data_processed/microbiome/data_asv_rra0.001_p0.02_th5000.csv")



######################################## old
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
tax <- read_delim("data/data_raw/data_microbiome/ASVs_taxonomy_new.tsv") %>% 
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
