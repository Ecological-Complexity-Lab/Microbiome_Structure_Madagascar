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
    #mutate(reads = 1) %>% 
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
  
  hosts <- dat %>% distinct(host_ID, grid, host_group) 
  
  nmi_obs <- aricode::NMI(hosts$grid, hosts$host_group)
  nmi_shuff <- vector(length = 1000)
  for(i in 1:1000) {
    # shuffling the grid attribute
    hosts_shuff <- hosts %>% 
      mutate(grid = sample(grid,nrow(hosts)))
    
    # calculating nmi
    nmi_shuff[i] <- aricode::NMI(hosts_shuff$grid, hosts_shuff$host_group)
  }
  
  # calculating p value
  p <- length(nmi_shuff[nmi_shuff>nmi_obs]) / length(nmi_shuff)
  
  # plotting
  if(figure){
    g <- as.data.frame(nmi_shuff) %>% 
      ggplot(aes(nmi_shuff)) + 
      geom_histogram() + 
      theme_classic() +
      geom_vline(xintercept = nmi_obs, linetype='dashed', color="red") +
      theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20)) +
      labs(x="NMI", y="No. of Shuffled Networks") +
      annotate(geom = "text", x=c((nmi_obs-0.01)), y=c(120,110), label=c(paste('NMI =',round(nmi_obs,3)), paste('p-value =',p)))
  } else {
    g <- NULL
  }
  
  nmi_summary <- tibble(nmi = nmi_obs,
                        p = p)
  
  return(list(nmi_summary, g))
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
    mutate(degree = i, type = "core", sig = ifelse(p <= 0.05, 1, 0))
  
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
    mutate(degree = j, type = "non-core", sig = ifelse(p <= 0.05, 1, 0))
  
  nmi_summary_noncore <- rbind(nmi_summary_noncore, nmi_mid)
}

#####
# plotting the results
nmi_summary <- rbind(nmi_summary_core, nmi_summary_noncore)

p1 <- nmi_summary %>% 
  ggplot(aes(x=degree, y=nmi, color = type)) + 
  geom_point() +
  geom_line() + 
  geom_hline(yintercept = nmi_observed[[1]]$nmi, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 0.35)) +
  scale_x_continuous(limits = c(0, 15)) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 18), plot.title = element_text(hjust = 0.5)) +
  labs(x="ASVs Degree", y="Normalized Mutual Information (NMI)")

return(list(p1))
}


# function for calculating similarity in modules between grids
fun_modules_similarity <- function(dat) {
  
  # matrix of grid similarity in modules
  grid_modules <- dat %>% 
    group_by(grid, host_group) %>% 
    summarise(host_n = n_distinct(host_ID)) %>% 
    spread(host_group, host_n, fill = 0) %>% 
    mutate(grid = as.character(grid)) %>%
    arrange(grid) %>% 
    column_to_rownames("grid") %>% 
    as.matrix() 
  
  # calculating the *similarity* between grids
  grid_modules_dist <- as.matrix(1-vegdist(sqrt(grid_modules), method = "bray"))
  
  # transforming to long format
  grid_mudules_dist_m <- grid_modules_dist
  grid_mudules_dist_m[lower.tri(grid_mudules_dist_m)] <- NA
  diag(grid_mudules_dist_m) <- NA
  grid_mudules_dist_m <- melt(grid_mudules_dist_m) %>% 
    filter(!is.na(value)) %>% 
    dplyr::rename(grid1 = Var2, grid2 = Var1, module_similarity = value) 
  
  return(grid_mudules_dist_m)
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
    ggplot(aes(x=n)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), strip.text.x = element_text(size=12)) +
    labs(x="Module Size", y="No. of Modules")
  
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

############################################################################

# ######################################################
# # hypotheses testing
# 

# 
# ##### model selection
# 
# # full model
# # adding the village as a random factor
# 
# 
# 
# library(nlme)
# full_model <- lme(module_similarity ~ grid_attr + grid_dist + sm_community, random=~1|village, data = final_data, method = "ML", na.action = na.fail)
# 
# #library(lme4)
# #full_model <- lme4::glmer(module_similarity ~ grid_attr + grid_dist + sm_community + (1|village), data = final_data, REML = F, na.action = na.fail, family = "beta")
# #library(glmmTMB)
# #full_model <- glmmTMB::glmmTMB(module_similarity ~ grid_attr + grid_dist + sm_community + (1|village), data = final_data, REML = F, na.action = na.fail, family = beta_family())
# 
# # checking VIF
# # as a rule of thumb, VIF< 10 for a variable is fine
# library(car)
# car::vif(full_model)
# 
# #library(fitdistrplus)
# #a=fitdist(final_data$module_similarity, "norm")
# #b=fitdist(final_data$module_similarity, "beta")
# #a$aic
# #b$aic
# #plot(b)
# #qqPlot(final_data$grid_attr)
# #shapiro.test(final_data$module_similarity)
# 




######################################################

# # number of modules as a function of asv degree
# n_module_grid <- modules_observed %>% 
#   group_by(host_group) %>% 
#   summarise(n_grid = n_distinct(grid))
# 
# n_asv <- modules_observed %>% 
#   #distinct(asv_ID, asv_degree, host_group) %>% 
#   left_join(n_module_grid, by = c("asv_group"="host_group")) %>% 
#   group_by(asv_ID, asv_degree) %>% summarise(n_grid = mean(n_grid)) %>% 
#   mutate(n = ifelse(asv_degree<5,"1-4", ifelse(asv_degree>=10,"10+", "5-9"))) %>% 
#   mutate(n = factor(n, levels = c("1-4","5-9","10+"))) 
# 
# p3 <- n_asv %>% 
#   ggplot(aes(x=n, y=n_grid, fill=n)) + 
#   geom_boxplot() + 
#   theme_bw() +
#   theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), legend.position = "none") +
#   labs(x="ASVs Degree", y="Modules' No. of Land Uses")
