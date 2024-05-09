# analyzing grid modularity using Infomap
# for different core microbes

library(tidyverse)
library(magrittr)
library(ggplot2)
library(vegan)
library(igraph)
library(reshape2)
library(infomapecology)
library(aricode)

rm(list=ls())

###########################################################################



##### functions #####################################################################

# function for finding modules
fun_modularity_analysis <- function(dat) {
  # input: dat = ASV data in long format
  # output: adding two columns to the input table: "asv_group" and "host_group"
  
  # transforming to matrix
  data_asv_mat <- dat %>% select(-village,-host_species,-grid,-season,-asv_degree,-total_reads) %>% 
    mutate(reads = 1) %>% 
    spread(asv_ID, reads, fill = 0) %>% 
    column_to_rownames("host_ID") %>% 
    as.matrix()
  
  # building the network
  network_object <- infomapecology::create_monolayer_object(data_asv_mat, directed = FALSE, bipartite = TRUE, group_names = c("ASV", "Host"))
  
  # modularity analysis
  infomap_object <- run_infomap_monolayer(network_object,
                                          infomap_executable='Infomap',
                                          flow_model = 'undirected',
                                          two_level = TRUE,
                                          silent=TRUE, trials=100, seed=123)
  
  # adding the modules classification to the data
  modules_asv <- infomap_object$modules %>% filter(node_group == "ASV")%>% dplyr::rename(asv_group = module_level1, asv_ID = node_name)
  modules_host <- infomap_object$modules %>% filter(node_group == "Host")%>% dplyr::rename(host_group = module_level1, host_ID = node_name) %>% mutate(host_ID = as.numeric(host_ID))
  
  dat %<>% left_join(modules_asv %>% dplyr::select(asv_ID,asv_group), by = "asv_ID") %>% 
    left_join(modules_host %>% dplyr::select(host_ID,host_group), by = "host_ID")
  
  return(dat)
}

# function for calculating NMI between the host's group and the host's grid and calculating p value for the value
fun_nmi_calc <- function(dat, figure) {
  # input: dat = ASV data in long format, figure = print the figure? (T/F)
  # output: data frame with observed NMI value and p value, if figure=T a figure of shuffled nmi distribution
  
  hosts <- dat %>% ungroup() %>% distinct(host_ID, grid, host_group)
  
  nmi_obs <- aricode::NMI(hosts$grid, hosts$host_group, "sum")
  nmi_shuff <- vector(length = 1000)
  for(i in 1:1000) {
    # shuffling the grid attribute
    hosts_shuff <- hosts %>% 
      mutate(grid = sample(grid,nrow(hosts)))
    
    # calculating nmi
    nmi_shuff[i] <- aricode::NMI(hosts_shuff$grid, hosts_shuff$host_group, "sum")
  }
  
  # calculating p value
  p <- length(nmi_shuff[nmi_shuff>=nmi_obs]) / length(nmi_shuff)
  
  # plotting
  if(figure){
    g <- as.data.frame(nmi_shuff) %>% 
      ggplot(aes(nmi_shuff)) + 
      geom_histogram(fill = "black") + 
      theme_bw() +
      geom_vline(xintercept = nmi_obs, linetype='dashed', color="red") +
      theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), panel.grid = element_blank()) +
      labs(x="Normalized Mutual Information (NMI)", y="Count") +
      annotate(geom = "text", x=c((nmi_obs-0.01)), y=c(110,100), label=c(paste('NMI =',round(nmi_obs,3)), paste('p-value =',p)))
  } else {
    g <- NULL
  }
  
  nmi_summary <- tibble(nmi = nmi_obs,
                        p = p)
  
  return(list(nmi_summary, nmi_shuff))
}


# calculating NMI for different values of core microbiome
# core microbiome = the number of hosts the asv infects (asv's degree)
fun_modularity_diff_core <- function(dat, nmi_observed, cor_seq) {
  
#####
# running for core microbiome
nmi_summary_core <- NULL
for (i in core_seq) {
  
  # filtering out ASVs with lower degree than the threshold
  data_asv_filtered_core <- dat %>% filter(asv_degree >= i)
  
  # modularity
  modules <- fun_modularity_analysis(data_asv_filtered_core)
  
  # NMI
  nmi_mid <- fun_nmi_calc(modules, FALSE)
  nmi_mid <- nmi_mid[[1]] %>% 
    mutate(degree = i, type = "core", sig = ifelse(p <= 0.05, 1, 0), n_module = length(unique(modules$host_group)))
  
  nmi_summary_core <- rbind(nmi_summary_core, nmi_mid)
}

#####
# running for non-core microbiome
nmi_summary_noncore <- NULL
for (j in core_seq) {
  
  # filtering out ASVs with lower degree than the threshold
  data_asv_filtered_noncore <- dat %>% filter(asv_degree <= j)
  
  # modularity
  modules <- fun_modularity_analysis(data_asv_filtered_noncore)
  
  # NMI
  nmi_mid <- fun_nmi_calc(modules, FALSE)
  nmi_mid <- nmi_mid[[1]] %>% 
    mutate(degree = j, type = "non-core", sig = ifelse(p <= 0.05, 1, 0), n_module = length(unique(modules$host_group)))
  
  nmi_summary_noncore <- rbind(nmi_summary_noncore, nmi_mid)
}

#####
# plotting the results
nmi_summary <- rbind(nmi_summary_core, nmi_summary_noncore)

p1 <- nmi_summary %>% 
  filter(type == "core") %>% 
  mutate(nmi = ifelse(is.na(nmi),0,nmi)) %>% 
  ggplot(aes(x=degree, y=nmi)) + 
  geom_hline(yintercept = nmi_observed[[1]]$nmi, linetype = "dashed") +
  geom_line(color = "black") +
  geom_point(aes(shape = as.factor(sig), size = n_module), color = "black") +
  scale_shape_manual(values = c(1, 16)) +
  scale_y_continuous(limits = c(0, 0.35)) +
  scale_x_continuous(limits = c(0, 20)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 15)) +
  labs(x="Minimum ASVs Degree", y="Normalized Mutual Information (NMI)")

return(nmi_summary)
}


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
fun_modules <- function(dat) {
  
  # how many host individuals in every grid
  n_host_grid <- dat %>% 
    distinct(host_ID, grid) %>% 
    count(grid) %>% 
    dplyr::rename(n_grid=n)
  
  g <- dat %>% 
    distinct(host_ID, grid, host_group) %>% 
    count(grid, host_group) %>% 
    left_join(n_host_grid, by="grid") %>% 
    mutate(n_rel = n/n_grid) %>% 
    ggplot(aes(x = host_group, y = grid, fill=n_rel)) +
    geom_tile(color='white') +
    theme_classic() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
    labs(x='Module ID', y='Land Use', color = "Host Abundance [%]")
  
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

fun_calc_betaNTI <- function(dat_mat, phylo_dist, asv_pool) {
  
  n_modules <- nrow(dat_mat)
  n_asv <- ncol(dat_mat)
  n_shuff <- 20
  
  # calculating observed MNTD
  mntd_obs <- as.matrix(picante::comdistnt(dat_mat, phylo_dist))
  
  mntd_shuff <- array(NA,dim = c(n_modules,n_modules,n_shuff))
  # loop for shuffling
  for (i in 1:n_shuff) {
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
  
  # transforming to long format
  betaNTI_mat2 <- betaNTI_mat
  betaNTI_mat2[upper.tri(betaNTI_mat2)] <- NA
  diag(betaNTI_mat2) <- NA
  betaNTI <- melt(betaNTI_mat2) %>% 
    filter(!(is.na(value))) %>% 
    rename(host1 = Var1, host2 = Var2, betaNTI = value)
  
  return(betaNTI)
}

############################################################################

