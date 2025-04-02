# functions for the scripts:
# main_analysis
# main_script

library(tidyverse)
library(magrittr)
library(ggplot2)
library(vegan)
library(igraph)
library(reshape2)
library(infomapecology)

rm(list=ls())



##### functions #####################################################################

##### 1
# function for finding modules
fun_modularity_analysis <- function(dat) {
  # input: dat = ASV data in long format
  # output: adding two columns to the input table: "asv_group" and "host_group"
  
  # transforming to matrix
  data_asv_mat <- dat %>% select(asv_ID, host_ID, reads) %>% 
    spread(asv_ID, reads, fill = 0) %>% 
    column_to_rownames("host_ID") %>% 
    as.matrix()
  
  # building the network
  network_object <- infomapecology::create_monolayer_network(data_asv_mat, directed = FALSE, bipartite = TRUE, group_names = c("ASV", "Host"))
  
  # modularity analysis
  infomap_object <- infomapecology::run_infomap_monolayer(network_object,
                                          infomap_executable='Infomap',
                                          flow_model = 'undirected',
                                          two_level = TRUE,
                                          silent=TRUE, trials=100, seed=123)
  
  # adding the modules classification to the data
  modules_asv <- infomap_object$modules %>% 
    filter(node_group == "ASV") %>% 
    dplyr::rename(asv_group = module_level1, asv_ID = node_name)
  
  modules_host <- infomap_object$modules %>% 
    filter(node_group == "Host") %>% 
    dplyr::rename(host_group = module_level1, host_ID = node_name) %>% 
    mutate(host_ID = as.numeric(host_ID))
  
  dat %<>% left_join(modules_asv %>% dplyr::select(asv_ID,asv_group), by = "asv_ID") %>% 
    left_join(modules_host %>% dplyr::select(host_ID,host_group), by = "host_ID")
  
  return(dat)
}


##### 2
# function for calculating similarity in modules between grids
# output: a matrix with sites as rows and modues as columns

fun_modules_similarity <- function(dat) {
  
  # matrix of grid similarity in modules
  grid_modules <- dat %>% 
    group_by(grid, village, host_group) %>% 
    summarise(n = n_distinct(host_ID)) %>% 
    spread(host_group, n, fill = 0) %>% 
    mutate(grid = as.character(grid)) %>%
    arrange(grid, village) %>% 
    unite("sample", c("village","grid"), remove = TRUE) %>% 
    column_to_rownames("sample") %>% 
    as.matrix() 
  
  return(list(grid_modules))
}


#####################################
#####################################
#####################################




# function for calculating similarity in modules between grids
fun_modules_similarity <- function(dat) {
  
  # how many host individuals in every grid
  n_host_grid <- dat %>% 
    distinct(host_ID, grid) %>% 
    count(grid) %>% 
    dplyr::rename(n_grid=n)
  
  # matrix of grid similarity in modules
  grid_modules <- dat %>% 
    group_by(grid, host_group) %>% 
    summarise(host_n = n_distinct(host_ID)) %>% 
    left_join(n_host_grid, by="grid") %>% 
    mutate(n = host_n/n_grid) %>% 
    select(-host_n, -n_grid) %>% 
    spread(host_group, n, fill = 0) %>% 
    mutate(grid = as.character(grid)) %>%
    arrange(grid) %>% 
    column_to_rownames("grid") %>% 
    as.matrix() 
  
  # calculating the *similarity* between grids
  grid_modules_dist <- as.matrix(1-vegdist(grid_modules, method = "bray"))
  
  # # transforming to long format
  # grid_mudules_dist_m <- grid_modules_dist
  # grid_mudules_dist_m[lower.tri(grid_mudules_dist_m)] <- NA
  # diag(grid_mudules_dist_m) <- NA
  # grid_mudules_dist_m <- melt(grid_mudules_dist_m) %>% 
  #   filter(!is.na(value)) %>% 
  #   dplyr::rename(grid1 = Var2, grid2 = Var1, module_similarity = value) 
  
  return(list(grid_modules_dist))
}



# function for plotting ASVs degree distribution
fun_asv_degree_distribution <- function(dat) {
  
  g <- dat %>% 
    distinct(asv_ID, asv_degree) %>% 
    ggplot(aes(x=asv_degree)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = c(10), linetype="dashed") +
    xlim(0,100)+
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
    labs(x="ASVs Degree", y="No. of ASVs")
  
  return(list(g))
}

# function for plotting modules
fun_modules <- function(dat, v) {
  
  # how many host individuals in every grid
  n_host_grid <- dat %>% 
    distinct(host_ID, village, grid) %>% 
    count(village, grid) %>% 
    dplyr::rename(n_grid=n) %>%
    mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth","agriculture","agroforest","flooded_rice","village"))) %>%
    arrange(grid) %>% 
    mutate(site = paste0("Site", row_number())) %>% 
    mutate(site = factor(site, levels = paste0("Site", row_number()))) 
  
  # color
  if (v == "Core"){
    asv_core_col1 <- "#d95f02"
    asv_core_col2 <- "#f8c49a"
  } else {
    asv_core_col1 <- "#034e7b"
    asv_core_col2 <- "#bbdefb"
  }
  
  # number of modules
  n_modules <- length(unique(dat$host_group))
  
  g <- dat %>% 
    distinct(host_ID, village, grid, host_group, asv_core) %>% 
    count(village, grid, host_group, asv_core) %>% 
    left_join(n_host_grid, by=c("village", "grid")) %>% 
    mutate(n_rel = n/n_grid) %>% 
    #mutate(host_group = as.factor(host_group)) %>% 
    unite("sample", c("village","grid"), remove=F) %>% 
    ggplot(aes(x = host_group, y = site, fill=n_rel)) +
    geom_tile(color='white') +
    theme_classic() +
    scale_x_continuous(breaks = seq(1,n_modules,length.out=5))+
    scale_fill_gradient(low = asv_core_col2, high = asv_core_col1, name = v) +
    theme(axis.text = element_text(size = 7, color = 'black'), title = element_text(size = 12), strip.text.x = element_text(size=11), 
          aspect.ratio = 0.8) +
    labs(title = paste(v), x='Module ID', y='', fill = "Host Relative Abundance")
 
  return(list(g))
}


# function for plotting modules size
fun_module_size <- function(dat) {
  
  g <- dat %>% 
    group_by(host_group) %>% 
    summarise(n = n_distinct(host_ID)) %>% 
    arrange(desc(n)) %>% 
    mutate(host_group = factor(host_group, levels = host_group)) %>% 
    ggplot(aes(x=factor(host_group), y=n)) +
    geom_bar(stat = "identity", width = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), axis.text.x=element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x="Modules", y="Module Size")
  
  return(list(g))
}


# function for plotting number of land uses per module
fun_module_grid <- function(dat) {
  
  n_landuse <- dat %>% 
    group_by(host_group) %>% 
    summarise(n = n_distinct(grid)) 
  
  g <- n_landuse %>%  
    ggplot(aes(x=n)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
    labs(x="No. of Land Uses", y="No. of Modules")
  
  return(list(g))
}


###########################################################################
# function to calculate beta-NTI

fun_calc_betaNTI <- function(dat_mat, phylo_dist, asv_pool, n_shuff) {
  
  n_modules <- nrow(dat_mat)
  n_asv <- ncol(dat_mat)
  seeds <- seq(1, n_shuff*100, length.out=n_shuff)
  
  # calculating observed MNTD
  mntd_obs <- as.matrix(picante::comdistnt(dat_mat, phylo_dist))
  
  mntd_shuff <- array(NA,dim = c(n_modules,n_modules,n_shuff))
  # loop for shuffling
  for (i in 1:n_shuff) {
    # randomly sampling ASVs from the ASV pool
    set.seed(seeds[i])
    shuff_asv_names <- sample(asv_pool$asv_ID, n_asv, replace = FALSE)
    #prob = asv_pool$p,
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
  
  # transforming to long format
  betaNTI_mat2 <- betaNTI_mat
  betaNTI_mat2[is.na(betaNTI_mat2)] <- 0
  betaNTI_mat2[upper.tri(betaNTI_mat2)] <- NA
  diag(betaNTI_mat2) <- NA
  betaNTI <- melt(betaNTI_mat2) %>% 
    filter(!(is.na(value))) %>% 
    dplyr::rename(host1 = Var1, host2 = Var2, betaNTI = value)
  
  return(betaNTI)
}
