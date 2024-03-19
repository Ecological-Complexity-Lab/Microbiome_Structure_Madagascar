# exploring microbes sharing 

library(tidyverse)
library(magrittr)
library(ggplot2)

rm(list=ls())


# reading microbiome data
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv")

# ASVs that have at least one link to host species other than rattus
asv_no_rattus <- data_asv %>% 
  filter(host_species != "Rattus rattus" & village == "Mandena") %>% 
  distinct(asv_ID) %>% 
  mutate(shared = 1)

dat <- data_asv %>% 
  filter(host_species == "Rattus rattus" & village == "Mandena")

# matrix for shared ASVs in the grid level
shared_mat <- dat %>% 
  filter(asv_ID %in% asv_no_rattus$asv_ID) %>% 
  group_by(grid, asv_ID) %>% 
  summarise(reads = mean(reads)) %>% 
  mutate(reads = 1) %>% 
  spread(asv_ID, reads, fill = 0) %>% 
  column_to_rownames("grid") %>% 
  dplyr::select(starts_with("ASV")) %>% 
  as.matrix()

shared_species <- as.matrix(1-vegdist(shared_mat, method = "jaccard"))

# long format
shared_species2 <- shared_species
shared_species2[upper.tri(shared_species2)] <- NA
diag(shared_species2) <- NA
grid_shared_similarity <- melt(shared_species2) %>% 
  rename(grid1=Var1, grid2=Var2, shared=value) %>% 
  filter(!(is.na(shared))) %>% 
  arrange(grid1, grid2)

final_data <- read_csv("data/data_processed/final_modularity_data.csv") %>% 
  filter(village == "Mandena") %>% 
  left_join(grid_shared_similarity, by=c("grid1","grid2"))

final_data %>% 
  ggplot(aes(y=shared, x=(1-sm_community))) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=T, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Small Mammals Similarity [Bray-Curtis]", y = "Shared ASVs similarity [Jaccard]")

##########################################
# matrix for all ASVs in the grid level
shared_all_mat <- dat %>% 
  filter(!(asv_ID %in% asv_no_rattus$asv_ID)) %>%
  group_by(grid, asv_ID) %>% 
  summarise(reads = mean(reads)) %>% 
  mutate(reads = 1) %>% 
  spread(asv_ID, reads, fill = 0) %>% 
  column_to_rownames("grid") %>% 
  dplyr::select(starts_with("ASV")) %>% 
  as.matrix()

no_shared_species <- as.matrix(1-vegdist(shared_all_mat, method = "jaccard"))

# long format
no_shared_species2 <- no_shared_species
no_shared_species2[upper.tri(no_shared_species2)] <- NA
diag(no_shared_species2) <- NA
grid_no_shared_similarity <- melt(no_shared_species2) %>% 
  rename(grid1=Var1, grid2=Var2, shared=value) %>% 
  filter(!(is.na(shared))) %>% 
  arrange(grid1, grid2)

final_data <- read_csv("data/data_processed/final_modularity_data.csv") %>% 
  filter(village == "Mandena") %>% 
  left_join(grid_shared_similarity, by=c("grid1","grid2"))

final_data %>% 
  ggplot(aes(y=shared, x=(1-sm_community))) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=T, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Small Mammals Similarity [Bray-Curtis]", y = "Shared ASVs similarity [Jaccard]")


##############################################
final_shared <- shared_species/no_shared_species

# long format
final_shared2 <- final_shared
final_shared2[upper.tri(final_shared2)] <- NA
diag(final_shared2) <- NA
grid_final_shared_similarity <- melt(final_shared2) %>% 
  rename(grid1=Var1, grid2=Var2, shared=value) %>% 
  filter(!(is.na(shared))) %>% 
  arrange(grid1, grid2)

final_data <- read_csv("data/data_processed/final_modularity_data.csv") %>% 
  filter(village == "Mandena") %>% 
  left_join(grid_final_shared_similarity, by=c("grid1","grid2"))

final_data %>% 
  ggplot(aes(y=shared, x=(1-sm_community))) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=T, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Small Mammals Similarity [Bray-Curtis]", y = "Shared ASVs similarity [Jaccard]")
