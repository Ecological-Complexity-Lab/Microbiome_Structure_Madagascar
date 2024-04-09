# analyzing phylogenetic distance between land uses

library(tidyverse)
library(magrittr)
library(ggplot2)
library(vegan)
library(ape)
library(picante)

rm(list=ls())

###########################################################################

# reading the phylogenetic tree
best_tree <- readRDS(file = "results/phylo_tree_0.01_2.rds")
phylo_tree <- best_tree$tree

# reading the microbiome data with modules
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000_with_modules.csv") %>% 
filter(host_species == "Rattus rattus" & grid!="village")


###########################################################################
# function to calculate beta-NTI

fun_calc_betaNTI <- function(dat_mat, phylo_dist, asv_pool) {
  
  n_modules <- nrow(dat_mat)
  n_asv <- ncol(dat_mat)
  
  #asv_pool2 <- colnames(dat_mat)
  
  # calculating observed MNTD
  mntd_obs <- as.matrix(picante::comdistnt(dat_mat, phylo_dist))
  
  mntd_shuff <- array(NA,dim = c(n_modules,n_modules,10))
  # loop for shuffling
  for (i in 1:10) {
    # randomly sampling ASVs from the ASV pool
    shuff_asv_names <- sample(asv_pool, n_asv, replace = FALSE)
    # changing the cols names
    dat_mat_shuff <- dat_mat
    colnames(dat_mat_shuff) <- shuff_asv_names
    # calculating MNTD
    mntd_shuff[,,i] <- as.matrix(picante::comdistnt(dat_mat_shuff, phylo_dist))
  }
  
  # taking the mean and sd of the shuffled matrices
  mntd_shuff_mean <- apply(mntd_shuff,c(1,2),mean)
  mntd_shuff_sd <- apply(mntd_shuff,c(1,2),sd)
  
  # calculating betaNTI
  betaNTI_mat <- (mntd_obs - mntd_shuff_mean) / mntd_shuff_sd
  betaNTI_mat2 <- betaNTI_mat
  betaNTI_mat2[upper.tri(betaNTI_mat2)] <- NA
  diag(betaNTI_mat2) <- NA
  betaNTI <- melt(betaNTI_mat2) %>% 
    filter(!(is.na(value))) %>% 
    rename(module1 = Var1, module2 = Var2, betaNTI = value)
    
  return(betaNTI)
}



###########################################################################
# main script

# ASVs phylogenetic distance
asv_distance <- ape::cophenetic.phylo(phylo_tree)

# filter villages
data_asv_village <- data_asv %>% 
  filter(village == "Mandena")

# ASVs pool
asv_pool <- unique(data_asv_village$asv_ID)

grids_names <- unique(data_asv_village$grid)
betaNTI_results <- NULL

# loop of grids
for (g in grids_names) {
  data_asv_grid <- data_asv_village %>% filter(grid == g)
  
  # transforming to a community matrix - rows = modules, cols = ASVs
  data_asv_mat <- data_asv_grid %>% 
    filter(host_group==asv_group) %>% 
    distinct(asv_ID, host_group) %>% 
    mutate(abundance = 1) %>% 
    spread(asv_ID, abundance, fill = 0) %>% 
    column_to_rownames("host_group") %>% 
    as.matrix()
  
  # calculating betaNTI
  betaNTI_results_grid <- fun_calc_betaNTI(data_asv_mat, asv_distance, asv_pool) %>% 
    mutate(grid = g)
  betaNTI_results <- rbind(betaNTI_results, betaNTI_results_grid)
}

# finding the small modules (occuring in only 1-2 grids)
n_grids_module <- data_asv_village %>% 
  group_by(host_group) %>% 
  summarise(n_grid = n_distinct(grid)) %>% 
  filter(n_grid <= 2)

# plotting
betaNTI_results %>% 
  ggplot(aes(x=grid, y=betaNTI, fill=grid)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text = element_text(size = 8, color = 'black'), title = element_text(size = 20), legend.position = "none") +
  labs(x="Land use", y="betaNTI")



###########################################################################

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
