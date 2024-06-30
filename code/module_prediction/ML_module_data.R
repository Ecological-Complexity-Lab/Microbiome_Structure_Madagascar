# script for making a table for ML 

library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)


rm(list=ls())


# villages: Andatsakala, Mandena, Sarahandrano
vil <- "Mandena"
group <- "Core"
final_table_three_villages <- NULL

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
  left_join(grid_similarity %>% select(grid1,grid2,grid_attr,sm_community), by=c("grid.x"="grid1", "grid.y"="grid2")) %>% 
  mutate(grid_attr = ifelse(grid.x==grid.y, 0, grid_attr), sm_community = ifelse(grid.x==grid.y, 0, sm_community)) %>% 
  mutate(elevation = abs(elevation.obs.x-elevation.obs.y)) %>% 
  mutate(module = ifelse(host_group.x==host_group.y, 1, 0)) %>% 
  select(-host_ID.x,-host_ID.y,-grid.x,-grid.y, -asv_core.x,-asv_core.y,-host_group.x,-host_group.y,-elevation.obs.x,-elevation.obs.y)



# saving table for three villages
final_table_three_villages <- rbind(final_table_three_villages, final_table)

# saving the final table as .csv
write_csv(final_table_three_villages, "data/data_processed/ML_module/ML_rattus_three_villages_core.csv")

