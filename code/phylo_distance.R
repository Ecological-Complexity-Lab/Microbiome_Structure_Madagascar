# analyzing phylogenetic distance between land uses

library(tidyverse)
library(magrittr)
library(ggplot2)
library(vegan)
library(ape)
library(phylocomr)

rm(list=ls())

###########################################################################

# reading the phylogenetic tree
best_tree <- readRDS(file = "results/phylo_tree_0.01_2.rds")
phylo_tree <- best_tree$tree
phylo_tree$tip.label <- tolower(phylo_tree$tip.label)
phylo_tree$node.label <- 1:length(phylo_tree$node.label)
phylo_tree_sing <- ape::collapse.singles(phylo_tree)

# reading the microbiome data
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv") %>% 
filter(host_species == "Rattus rattus" & grid!="village" & village == "Mandena")

# preparing the input table for the beta-NTI function
data_betaNTI <- data_asv %>% 
  group_by(grid, asv_ID) %>% 
  summarise(reads = mean(reads)) %>% 
  mutate(abundance = 1) %>% 
  #mutate(asv_ID = gsub('_','', asv_ID)) %>% 
  mutate(asv_ID = tolower(asv_ID)) %>%
  dplyr::select(grid, abundance, asv_ID) %>% 
  as.data.frame()
  
options(scipen = 999)
# calculating betaNTI between grids
betaNTI_resutlts <- phylocomr::ph_comdistnt(data_betaNTI, phylo_tree_sing, rand_test = TRUE, null_model = 1, randomizations = 999, abundance = FALSE)
