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
  
  mntd_shuff <- array(NA,dim = c(n_modules,n_modules,20))
  # loop for shuffling
  for (i in 1:20) {
    # randomly sampling ASVs from the ASV pool
    shuff_asv_names <- sample(asv_pool$asv_ID, n_asv, prob = asv_pool$p, replace = FALSE)
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
asv_pool <- data_asv_village %>% 
  distinct(asv_ID, asv_degree) %>% 
  mutate(p = asv_degree/length(unique(data_asv_village$asv_ID)))

#asv_pool <- unique(data_asv_village$asv_ID)

grids_names <- unique(data_asv_village$grid)
betaNTI_results <- NULL

# loop of grids
for (g in grids_names) {
  data_asv_grid <- data_asv_village %>% filter(grid == g)
  
  # transforming to a community matrix - rows = modules, cols = ASVs
  data_asv_mat <- data_asv_grid %>% 
    filter(host_group==asv_group) %>% 
    filter(asv_degree <= 2) %>%
    distinct(asv_ID, host_group) %>% 
    mutate(abundance = 1) %>% 
    spread(asv_ID, abundance, fill = 0) %>% 
    column_to_rownames("host_group") %>% 
    as.matrix()
  
  ###
  # transforming to a community matrix - rows = individuals, cols = ASVs
  data_asv_mat2 <- data_asv_grid %>% 
    mutate(abundance = 1) %>% 
    #filter(asv_degree > 5) %>% 
    distinct(host_ID, asv_ID, abundance) %>% 
    spread(asv_ID, abundance, fill = 0) %>% 
    column_to_rownames("host_ID") %>% 
    as.matrix()
  
  # calculating betaNTI
  betaNTI_results_grid <- fun_calc_betaNTI(data_asv_mat2, asv_distance, asv_pool, FALSE) %>% 
    mutate(grid = g)
  betaNTI_results <- rbind(betaNTI_results, betaNTI_results_grid)
  
}

# finding the small modules (occuring in only 1-2 grids)
n_grids_module <- data_asv_village %>% 
  group_by(host_group) %>% 
  summarise(n_grid = n_distinct(grid)) 

# plotting
a=betaNTI_results %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group), by=c("module1"="host_ID")) %>% rename(host_group1=host_group) %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group), by=c("module2"="host_ID")) %>% rename(host_group2=host_group) %>%
  #filter(host_group1 %in% n_grids_module$n_grid & host_group2 %in% n_grids_module$n_grid) %>% 
  filter(host_group1 != host_group2) %>% 
  group_by(grid) %>% 
  summarise(mean = mean(betaNTI))
  group_by(grid) %>% 
  summarise(mean = mean(mean))
  ggplot(aes(x=grid, y=betaNTI, fill=grid)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text = element_text(size = 8, color = 'black'), title = element_text(size = 20), legend.position = "none") +
  labs(x="Land use", y="betaNTI")

b=a %>% filter(grid=="flooded_rice") %>% filter(host_group1 != host_group2)
t.test(b$mean, mu=0)


################
# individuals across grids
data_betaNTI <- data_asv_village %>% 
  distinct(host_ID, asv_ID) %>% 
  mutate(reads = 1) %>% 
  spread(asv_ID, reads, fill = 0) %>% 
  column_to_rownames("host_ID") %>% 
  as.matrix()

r <- fun_calc_betaNTI(data_betaNTI, asv_distance, asv_pool)

# the same module
same_module <- r %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module1"="host_ID")) %>% rename(host_group1=host_group, grid1=grid) %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module2"="host_ID")) %>% rename(host_group2=host_group, grid2=grid) %>%
  filter(host_group1 != host_group2) %>% mutate(sig = ifelse(betaNTI>2,1,0))
  group_by(host_group1) %>% 
  summarise(mean = mean(betaNTI)) %>% 
    left_join(n_grids_module, by=c("host_group1"="host_group"))

plot(same_module$mean~same_module$n_grid)
t.test(same_module$mean, mu=0)

raupc <- raupcrick(data_betaNTI)
raupc2 <- as.matrix(raupc) 
raupc2[upper.tri(raupc2)] <- NA
diag(raupc2) <- NA
same_module_raup <- melt(raupc2) %>% 
  filter(!(is.na(value)))%>% 
  rename(module1 = Var1, module2 = Var2, betaNTI = value) %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module1"="host_ID")) %>% rename(host_group1=host_group, grid1=grid) %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module2"="host_ID")) %>% rename(host_group2=host_group, grid2=grid) %>%
  filter(host_group1 == host_group2) %>% 
  group_by(host_group1) %>% 
  summarise(mean = mean(betaNTI))

################
# population level
n_grid <- data_asv_village %>% group_by(grid) %>% summarise(n_host = n_distinct(host_ID))

data_betaNTI <- data_asv_village %>% 
  group_by(grid, asv_ID) %>% 
  summarise(n = n_distinct(host_ID)) %>% 
  left_join(n_grid, by="grid") %>% 
  mutate(n_rel = n/n_host) %>% 
  select(-n,-n_host) %>% 
  spread(asv_ID, n_rel, fill = 0) %>% 
  column_to_rownames("grid") %>% 
  as.matrix()

betaNTI_pop <- fun_calc_betaNTI(data_betaNTI, asv_distance, asv_pool, TRUE)

###
across_grids <- r %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module1"="host_ID")) %>% rename(host_group1=host_group, grid1=grid) %>% 
  left_join(data_asv_village %>% distinct(host_ID, host_group, grid), by=c("module2"="host_ID")) %>% rename(host_group2=host_group, grid2=grid) %>%
  filter(grid1!=grid2) %>% 
  mutate(grid1_org = grid1) %>% 
  rowwise() %>%
  mutate(grid1 = sort(c(grid1, grid2))[1], grid2 = sort(c(grid1_org, grid2))[2]) %>% 
  group_by(grid1, grid2) %>% 
  summarise(mean = mean(betaNTI))

# reading the final modularity analysis dat
final_data <- read_csv("data/data_processed/final_modularity_data.csv") %>% 
  filter(village=="Mandena") %>% 
  left_join(across_grids, by=c("grid1"="grid2", "grid2"="grid1"))

final_data %>% 
  ggplot(aes(y=mean, x=module_similarity)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=F, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Small Mammals Disimilarity [Bray-Curtis]", y = "betaNTI")

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
  
# calculating betaNTI between grids
betaNTI_resutlts <- phylocomr::ph_comdistnt(data_betaNTI, phylo_tree_sing, rand_test = TRUE, null_model = 1, randomizations = 999, abundance = FALSE)
