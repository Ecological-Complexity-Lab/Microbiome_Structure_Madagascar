# script for making a table for ML 

library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)


rm(list=ls())


vil <- "Mandena"
group <- "Core"

#####################################################
# reading the data

# small mammals (hosts)
data_mammals_full <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv") %>% 
 mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) 

 data_mammals <- data_mammals_full %>% 
   select(host_ID, elevation.obs, sex, mass, age_repro) %>% 
  mutate(sex = as.factor(sex), age_repro = as.factor(age_repro)) 

# host modules
data_host <- read_csv("results/modules_table_three_villages.csv")

host_richness <- data_host %>% 
  group_by(host_ID, asv_core) %>% 
  summarise(richness = n_distinct(asv_ID))

data_host_filtered <- data_host %>% 
  filter(asv_core == group) %>% 
  distinct(host_ID, grid, season, host_group, asv_core) %>% 
  left_join(data_mammals, by="host_ID") %>% 
  left_join(host_richness, by=c("asv_core","host_ID")) %>% 
  mutate(season = factor(season))


#####################################################
# vegetation attributes and small mammals community *dis-similarity*

grid_similarity <- read_csv("data/data_processed/village_summary.csv") %>% 
  filter(village == vil) %>% 
  select(grid1, grid2, grid_attr, sm_community)

grid_similarity2 <- grid_similarity %>% rename(grid1=grid2, grid2=grid1)
grid_similarity <- rbind(grid_similarity, grid_similarity2)


#####################################################
# distance between hosts

data_host_dist <- data_host_filtered %>% 
  select(host_ID) %>% 
  left_join(data_mammals_full %>% select(host_ID, longitude, latitude), by="host_ID")

# calculating distance
library(geosphere)
host_dist <- geosphere::distm(data_host_dist[2:3] , fun = distHaversine)
rownames(host_dist) <- data_host_dist$host_ID
colnames(host_dist) <- data_host_dist$host_ID

# transforming to long format
host_dist2 <- host_dist
host_dist2[upper.tri(host_dist2)] <- NA
diag(host_dist2) <- NA
host_distance_m <- melt(host_dist2) %>% 
  rename(host_ID.x=Var1, host_ID.y=Var2, distance=value) %>% 
  filter(!(is.na(distance))) 

#####################################################
# making the full table

final_table <- host_distance_m %>% 
  left_join(data_host_filtered, by=c("host_ID.x"="host_ID")) %>% 
  left_join(data_host_filtered, by=c("host_ID.y"="host_ID")) %>% 
  left_join(grid_similarity %>% select(grid1,grid2,grid_attr,sm_community), by=c("grid.x"="grid1", "grid.y"="grid2"))








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

