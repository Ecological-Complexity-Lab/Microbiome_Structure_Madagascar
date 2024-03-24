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

# microbiome
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv")
data_asv_filtered <- data_asv %>% 
  filter(host_species == "Rattus rattus" & village == "Mandena" & grid != "village") %>% 
  select(host_ID, asv_ID, grid, season) %>% 
  left_join(data_mammals, by="host_ID") %>% 
  mutate(season = factor(season))

# asv_degree
asv_degree <- data_asv_filtered %>% 
  group_by(asv_ID) %>% 
  summarise(n = n_distinct(host_ID)) 
data_asv_filtered %<>% left_join(asv_degree, by="asv_ID") %>% 
  filter(n>1) %>% 
  select(-n)

#####################################################
# ASVs taxonomy

asv_taxa <- read_delim("data/data_raw/data_microbiome/ASVs_all_merged_taxonomy.tsv") %>%
  select(ASV, Phylum, Class)


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
# calculate network attributes

library(bipartite)

# transforming to matrix
data_asv_mat <- data_asv_filtered %>% select(host_ID, asv_ID) %>% 
  mutate(reads = 1) %>% 
  spread(asv_ID, reads, fill = 0) %>% 
  column_to_rownames("host_ID") %>% 
  as.matrix()

# centrality
node_centrality <- bipartite::specieslevel(data_asv_mat, index = c("degree", "betweenness", "closeness"))
asv_centrality <- node_centrality$`higher level` %>% 
  select(degree, betweenness, closeness) %>% 
  rename(asv_degree=degree, asv_betweenness=betweenness, asv_closeness=closeness) %>% 
  rownames_to_column("asv_ID")
host_centrality <- node_centrality$`lower level` %>% 
  select(degree, betweenness, closeness) %>% 
  rename(host_degree=degree, host_betweenness=betweenness, host_closeness=closeness) %>%
  rownames_to_column("host_ID") %>% mutate(host_ID = as.double(host_ID))

# modularity
# building the network
network_object <- infomapecology::create_monolayer_object(data_asv_mat, directed = FALSE, bipartite = TRUE, group_names = c("ASV", "Host"))

# modularity analysis
infomap_object <- run_infomap_monolayer(network_object,
                                        infomap_executable='Infomap',
                                        flow_model = 'undirected',
                                        two_level = TRUE,
                                        silent=TRUE, trials=100, seed=123)

# asv modules
modules_asv <- infomap_object$modules %>% 
  filter(node_group == "ASV") %>% 
  dplyr::rename(asv_module = module_level1, asv_ID = node_name) %>% 
  select(asv_ID, asv_module)
# host modules
modules_host <- infomap_object$modules %>% 
  filter(node_group == "Host") %>% 
  dplyr::rename(host_module = module_level1, host_ID = node_name) %>% 
  mutate(host_ID = as.numeric(host_ID)) %>% 
  select(host_ID, host_module)

# shortest path
short_path <- igraph::distances(network_object$igraph_object) %>% 
  as.data.frame() %>% 
  select(starts_with("ASV")) %>% 
  rownames_to_column("host_ID") %>% 
  mutate(host_ID = as.double(host_ID)) %>% 
  filter(host_ID %in% modules_host$host_ID) %>% 
  gather("asv_ID", "shortest_path",starts_with("ASV"))


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
write_csv(data_asv_mat, "data/data_processed/microbiome/ML_microgale_mandena.csv")

