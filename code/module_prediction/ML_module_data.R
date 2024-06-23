# script for making a table for ML 

library(tidyverse)
library(magrittr)


rm(list=ls())

#####################################################
# reading the data

# small mammals (hosts)
data_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")
data_mammals %<>% mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>% 
  select(host_ID, elevation.obs, sex, mass, age_repro) %>% 
  mutate(sex = as.factor(sex), age_repro = as.factor(age_repro)) 

# host modules
data_host <- read_csv("results/modules_table_three_villages.csv")

host_richness <- data_host %>% 
  group_by(host_ID, asv_core) %>% 
  summarise(richness = n_distinct(asv_ID))

data_host_filtered <- data_host %>% 
  filter(asv_core == "Core") %>% 
  distinct(host_ID, grid, season, host_group, asv_core) %>% 
  left_join(data_mammals, by="host_ID") %>% 
  left_join(host_richness, by=c("asv_core","host_ID")) %>% 
  mutate(season = factor(season))





#####################################################
# PCA for grid attributes

grid_attributes <- read_csv("data/data_processed/village_attributes.csv") %>% 
  rename(grid=grid_name) %>% 
  filter(village == "Mandena")

grid_attr_pca <- prcomp(grid_attributes[3:10], scale=TRUE)
# taking the value of pca1
attr_pca1 <- grid_attr_pca$x[,1]
grid_attributes %<>% mutate(pca_grid_attr = attr_pca1) %>% 
  select(grid, pca_grid_attr)

#####################################################
# PCA for small mammals community
data_mammals_full <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")

data_mammals_full_mat <- data_mammals_full %>% 
  filter(grepl("TMR", animal_id)) %>% 
  rename(host_species = field_identification, grid = habitat_type) %>% 
  filter(host_species != "Rattus rattus", grid != "village") %>% 
  select(host_species, grid) %>% 
  count(grid, host_species) %>% 
  spread(host_species, n, fill = 0) 

grid_sm_pca <- prcomp(data_mammals_full_mat[-1], scale=TRUE)
# taking the value of pca1
sm_pca1 <- grid_sm_pca$x[,1]

data_mammals_full_mat %<>% mutate(pca_grid_sm = sm_pca1) %>% 
  select(grid, pca_grid_sm)




#####################################################
# making the full table

final_table <- data_asv_filtered %>% 
  mutate(link = 1) %>% 
  spread(asv_ID, link, fill = 0)  %>% 
  gather("asv_ID","Link", starts_with("ASV")) %>% 
  left_join(grid_attributes, by="grid") %>% 
  left_join(data_mammals_full_mat, by="grid") %>% 
  left_join(asv_taxa, by=c("asv_ID"="ASV")) %>% 
  left_join(host_centrality, by="host_ID") %>% 
  left_join(asv_centrality, by="asv_ID") %>% 
  left_join(modules_host, by="host_ID") %>% 
  left_join(modules_asv, by="asv_ID") %>% 
  left_join(short_path, by=c("host_ID", "asv_ID")) %>% 
  mutate(pref_attach = host_degree*asv_degree) %>% 
  select(host_ID, asv_ID, Link, grid, pca_grid_attr, pca_grid_sm, season,
         elevation.obs, mass, sex, age_repro, Phylum, Class,
         host_degree, host_betweenness, host_closeness, host_module,
         asv_degree, asv_betweenness, asv_closeness, asv_module,
         shortest_path, pref_attach)





# saving the final table as .csv
write_csv(final_table, "ML_rattus_mandena_include_degree1.csv")

