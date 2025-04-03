# functions for the scripts:
# main_analysis_continuous
# main_analysis_threshold

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
  network_object <- create_monolayer_network(data_asv_mat, directed = FALSE, bipartite = TRUE, group_names = c("ASV", "Host"))
  
  # modularity analysis
  infomap_object <- run_infomap_monolayer(network_object,
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


##### 3
# function for plotting ASVs degree distribution
fun_asv_degree_distribution <- function(dat) {
  
  g <- dat %>% 
    distinct(asv_ID, asv_degree) %>% 
    ggplot(aes(x=asv_degree)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = c(10), linetype="dashed") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = 'black'), 
          title = element_text(size = 14), 
          strip.text.x = element_text(size=12)
    ) +
    labs(x="ASVs Degree", y="No. of ASVs")
  
  return(list(g))
}


##### 4
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
  
  # colors
  if (v == "Core"){
    asv_core_col1 <- "#d95f02"
    asv_core_col2 <- "#f8c49a"
  } else {
    asv_core_col1 <- "#034e7b"
    asv_core_col2 <- "#bbdefb"
  }
  
  # number of modules
  n_modules <- length(unique(dat$host_group))
  
  # plotting
  g <- dat %>% 
    distinct(host_ID, village, grid, host_group, asv_core) %>% 
    count(village, grid, host_group, asv_core) %>% 
    left_join(n_host_grid, by=c("village", "grid")) %>% 
    mutate(n_rel = n/n_grid) %>% 
    unite("sample", c("village","grid"), remove=F) %>% 
    ggplot(aes(x = host_group, y = site, fill=n_rel)) +
    geom_tile(color='white') +
    theme_classic() +
    scale_x_continuous(breaks = seq(1,n_modules,length.out=5))+
    scale_fill_gradient(low = asv_core_col2, high = asv_core_col1, name = v) +
    theme(axis.text = element_text(size = 7, color = 'black'), 
          title = element_text(size = 12), 
          strip.text.x = element_text(size=11), 
          aspect.ratio = 0.8
    ) +
    labs(title = paste(v), x='Module ID', y='', fill = "Host Relative Abundance")
  
  return(list(g))
}


##### 5
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
    theme(axis.text = element_text(size = 14, color = 'black'), 
          title = element_text(size = 20), 
          axis.text.x=element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    labs(x="Modules", y="Module Size")
  
  return(list(g))
}


##### 6
# function for plotting number of land uses per module
fun_module_grid <- function(dat) {
  
  n_landuse <- dat %>% 
    group_by(host_group) %>% 
    summarise(n = n_distinct(grid)) 
  
  g <- n_landuse %>%  
    ggplot(aes(x=n)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, color = 'black'), 
          title = element_text(size = 20), 
          strip.text.x = element_text(size=12)
    ) +
    labs(x="No. of Land Uses", y="No. of Modules")
  
  return(list(g))
}
