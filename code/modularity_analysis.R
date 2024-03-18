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
library(cowplot)

rm(list=ls())

############################################################################

# reading microbiome data
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv") %>% 
  filter(host_species == "Rattus rattus")



##### functions #####################################################################

# function for finding modules
fun_modularity_analysis <- function(dat) {
  # input: dat = ASV data in long format
  # output: adding two columns to the input table: "asv_group" and "host_group"
  
  # transforming to matrix
  data_asv_mat <- dat %>% select(-asv_degree) %>% 
    mutate(reads = 1) %>% 
    spread(asv_ID, reads, fill = 0) %>% 
    column_to_rownames("host_ID") %>% 
    dplyr::select(starts_with("ASV")) %>% 
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
      theme_bw() +
      geom_vline(xintercept = nmi_obs, linetype='dashed', color="red") +
      theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20)) +
      labs(x="NMI", y="No. of Shuffled Networks")
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
#nmi_summary <- rbind(nmi_summary_core, nmi_summary_noncore)

p1 <- nmi_summary_core %>% 
  ggplot(aes(x=degree, y=nmi)) + 
  geom_point(color = "blue") +
  geom_line(color = "blue") + 
  geom_hline(yintercept = nmi_observed[[1]]$nmi, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 0.35)) +
  scale_x_continuous(limits = c(0, 15)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 14), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  labs(title = "Core Microbes", x="Minimum ASVs Degree", y="Normalized Mutual Information (NMI)")

p2 <- nmi_summary_noncore %>% 
  ggplot(aes(x=degree, y=nmi)) + 
  geom_point(color = "red") +
  geom_line(color="red") + 
  geom_hline(yintercept = nmi_observed[[1]]$nmi, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 0.35)) +
  scale_x_continuous(limits = c(0, 15)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'), title = element_text(size = 14), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  labs(title = "Non-Core Microbes", x="Maximum ASVs Degree", y="Normalized Mutual Information (NMI)")

p3 <- cowplot::plot_grid(p1,p2)
final_figs <- list(p3)

return(final_figs)
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




############################################################################
# main script


# setting thresholds for core
core_seq <- seq(1:5)
village_names <- unique(data_asv$village)
nmi_observed_three_villages <- NULL
nmi_diff_core_three_villages <- NULL
modules_similarity_three_villages <- NULL

# for loop for three villages
for (v in village_names) {
  
  data_asv_village <- data_asv %>% 
    filter(village == v)

# calculating ASVs degree 
asv_degree <- data_asv_village %>% 
  count(asv_ID) %>% 
  rename(asv_degree = n)
data_asv_village %<>% left_join(asv_degree, by="asv_ID")

##### observed network
# finding modules
modules_observed <- fun_modularity_analysis(data_asv_village)

# calculating similarity in modules between grids
modules_similarity <- fun_modules_similarity(modules_observed) %>% 
  mutate(village = v)
# saving results in one table
modules_similarity_three_villages <- rbind(modules_similarity_three_villages, modules_similarity)

# calculating NMI
nmi_observed <- fun_nmi_calc(modules_observed, TRUE)
# saving results in one list
nmi_observed_three_villages <- append(nmi_observed_three_villages, nmi_observed)

##### NMI for different values of core microbiome
# calling the function
nmi_diff_core <- fun_modularity_diff_core(data_asv_village, nmi_observed, core_seq)
# saving results in one list
nmi_diff_core_three_villages <- append(nmi_diff_core_three_villages, nmi_diff_core)
}


######################################################
# hypotheses testing

# reading grid similarity results
grids_similarity_attr <- read_csv("data/data_processed/village_summary.csv")

# combining the variables
final_data <- modules_similarity_three_villages %>% 
  left_join(grids_similarity_attr, by=c("village","grid1","grid2")) %>% 
  mutate(grid_dist = sqrt(grid_dist)) %>% 
  filter(grid1!="village") # removing the village grid

# saving the results
write_csv(final_data, "data/data_processed/final_modularity_data.csv")


final_data <- read_csv("data/data_processed/final_modularity_data.csv")
##### plotting the regression
# grid attributes
final_data %>% 
  ggplot(aes(y=module_similarity, x=grid_attr, color=village)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=F, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Grids Attributes Disimilarity [Bray-Curtis]", y = "Modules Similarity [Bray-Curtis]")

# grid distance
final_data %>% 
  ggplot(aes(y=module_similarity, x=grid_dist, color=village)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=F, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Distance Between Grids [Log(m)]", y = "Modules Similarity [Bray-Curtis]")

# small mammals similarity
final_data %>% 
  ggplot(aes(y=module_similarity, x=sm_community, color=village)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "glm", se=F, method.args = list(family = "gaussian")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'), title = element_text(size = 14), strip.text.x = element_text(size=12)) +
  labs(x = "Small Mammals Disimilarity [Bray-Curtis]", y = "Modules Similarity [Bray-Curtis]")


##### model selection

# full model
# adding the village as a random factor

#library(lme4)
#full_model <- lme4::glmer(module_similarity ~ grid_attr + grid_dist + sm_community + (1|village), data = final_data, REML = F, na.action = na.fail, family = "beta")

library(nlme)
full_model <- lme(module_similarity ~ grid_attr + grid_dist + sm_community, random=~1|village, data = final_data, method = "ML", na.action = na.fail)

#library(glmmTMB)
#full_model <- glmmTMB::glmmTMB(module_similarity ~ grid_attr + grid_dist + sm_community + (1|village), data = final_data, REML = F, na.action = na.fail, family = beta_family())

# checking VIF
# as a rule of thumb, VIF< 10 for a variable is fine
library(car)
car::vif(full_model)

#library(fitdistrplus)
#a=fitdist(final_data$module_similarity, "norm")
#b=fitdist(final_data$module_similarity, "beta")
#a$aic
#b$aic
#plot(b)
qqPlot(final_data$grid_attr)
shapiro.test(final_data$module_similarity)

# AIC 
library(MuMIn)
dredge_modules_similarity <- MuMIn::dredge(full_model)
# The best models: delta <= 10
results_modules_similarity <- subset(dredge_modules_similarity, delta <= 10 | df == 3 | df == max(df), recalc.weights = FALSE)
row.names(results_modules_similarity) <- c(1:(nrow(results_modules_similarity)))

# variables importance
imp <- as.data.frame(MuMIn::sw(dredge_modules_similarity))
var_names <- rownames(imp)
imp_values <- as.vector(imp[[1]])
# Plotting
barplot(imp_values, names.arg=var_names, ylab = "Importance", ylim = c(0,1))





######################################################

# number of modules as a function of asv degree
n_module_grid <- modules_observed %>% 
  group_by(host_group) %>% 
  summarise(n_grid = n_distinct(grid))

n_asv <- modules_observed %>% 
  #distinct(asv_ID, asv_degree, host_group) %>% 
  left_join(n_module_grid, by = c("asv_group"="host_group")) %>% 
  group_by(asv_ID, asv_degree) %>% summarise(n_grid = mean(n_grid)) %>% 
  mutate(n = ifelse(asv_degree<5,"1-4", ifelse(asv_degree>=10,"10+", "5-9"))) %>% 
  mutate(n = factor(n, levels = c("1-4","5-9","10+"))) 

p3 <- n_asv %>% 
  ggplot(aes(x=n, y=n_grid, fill=n)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), legend.position = "none") +
  labs(x="ASVs Degree", y="Modules' No. of Land Uses")
