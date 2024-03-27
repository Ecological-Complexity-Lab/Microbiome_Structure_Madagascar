# script for analyze grids attributes
# *taken from Kayla script*

library(tidyverse)
library(magrittr)
library(vegan)
library(geosphere)
library(reshape2)
library(psych)

rm(list=ls())

############################################################################
# functions

##### grid attributes #######################################################################

fun_grid_attributes <- function(dat) {
  

habitat <- dat %>% 
  mutate(grid_name = case_when(grid_name=="Primary forest" ~ "semi-intact_forest",
                               grid_name=="Degraded forest" ~ "semi-intact_forest",
                               grid_name=="Secondary forest" ~ "secondary_forest",
                               grid_name=="Savoka" ~ "brushy_regrowth",
                               grid_name=="Sugar cane" ~ "agriculture",
                               grid_name=="Mixed agricultural fields" ~ "agriculture",
                               grid_name=="Vanilla, coffee, banana fields" ~ "agriculture",
                               grid_name=="Rice field" ~ "flooded_rice",
                               grid_name=="Vanilla" ~ "agroforest",
                               grid_name=="Village" ~ "village")) %>% 
  select(village, grid_name, plot, square, dead_logs, tree_dbh, tree_height, #omiting crown distance
         liana_number, herbaceous_height, herbaceous_cover, canopy_cover) #omitting observations and season
## deal with multitrunk trees by taking sqrt of sum of square DBH of each stem
habitat <- habitat %>% 
  separate(tree_dbh, c("dbh1", "dbh2", "dbh3"), sep = ",") %>% #using unique(habitat$tree_dbh) see max of 3 stems, give each their own column
  mutate(across(c(dbh1, dbh2, dbh3), as.numeric), #make numeric
         across(c(dbh1, dbh2, dbh3), ~ifelse(is.na(.x), 0, .x)), #make na into 0
         tree_dbh = dbh1^2+dbh2^2+dbh3^2, #sum sq dbhs
         tree_dbh = ifelse(tree_dbh == 0, NA, sqrt(tree_dbh))) %>% #sqrt of sum, if 0 make NA (no trees)
  select(-dbh1, -dbh2, -dbh3)
## deal with the herbaceous_height "x to x"
habitat <- habitat %>% 
  mutate(herbaceous_height = case_when(herbaceous_height == "0.4 to 1.5" ~ mean(c(0.4, 1.5)),
                                       herbaceous_height == "1.5 to 2" ~ mean(c(1.5, 2)),
                                       herbaceous_height == "0.3 to 1" ~ mean(c(0.3, 1)), 
                                       herbaceous_height == "0,3 to 1" ~ mean(c(0.3, 1)), 
                                       herbaceous_height == "0 to 1" ~ mean(c(0, 1)),
                                       herbaceous_height == "0 to 2" ~ mean(c(0, 2)), 
                                       TRUE ~ as.numeric(herbaceous_height)))
# add pitfall column
habitat <- habitat %>% 
  mutate(pitfall = ifelse(str_detect(plot, "Pitfall"), "Pitfall", "Grid"),
         pitfall_line = case_when(str_detect(plot, "Pitfall 1") ~ 1,
                                  str_detect(plot, "Pitfall 2") ~ 2,
                                  str_detect(plot, "Pitfall 3") ~ 3,
                                  str_detect(plot, "Pitfall 4") ~ 4,
                                  str_detect(plot, "Pitfall 5") ~ 5,
                                  str_detect(plot, "Pitfall 6") ~ 6,
                                  str_detect(plot, "Pitfall 7") ~ 7,
                                  str_detect(plot, "Pitfall 8") ~ 8,
                                  str_detect(plot, "Pitfall 9") ~ 9,
                                  str_detect(plot, "Pitfall 10") ~ 10,
                                  str_detect(plot, "Pitfall 11") ~ 11,
                                  str_detect(plot, "Pitfall 12") ~ 12,
                                  TRUE ~ 0))

## make a trees df
trees <- habitat %>% 
  select(village, grid_name, pitfall, plot, square, tree_height, tree_dbh) %>% 
  group_by(village, grid_name, plot, square) %>%
  summarise(tree_height = mean(tree_height, na.rm=T),
          tree_dbh = mean(tree_dbh, na.rm=T))

## create a new df individual tree measurements
squares <- habitat %>% 
  mutate(tree = ifelse(is.na(tree_dbh), 0, 1)) %>% ## count the number of trees in each square
  group_by(village, grid_name, plot, square) %>% 
  mutate(n_trees = sum(tree)) %>% 
  ungroup() %>% 
  select(-tree_height, -tree_dbh, -tree) %>% # remove tree info from habitat
  distinct() %>% 
  left_join(trees, by=c("village","grid_name","plot","square"))

# summarizing for plots
plot_sum <- squares %>% 
  mutate(across(c(n_trees, tree_height, tree_dbh, dead_logs, 
                  liana_number, herbaceous_height, 
                  herbaceous_cover, canopy_cover))) %>% 
  group_by(pitfall, pitfall_line, village, grid_name, plot) %>% 
  summarise(n_trees = sum(n_trees, na.rm=T),
            tree_height = mean(tree_height, na.rm=T),
            tree_dbh = mean(tree_dbh, na.rm=T),
            n_logs = sum(dead_logs, na.rm=T),
            n_liana = sum(liana_number, na.rm=T),
            m_herb_ht = mean(herbaceous_height, na.rm=T),
            m_herb_cv = mean(herbaceous_cover, na.rm = T),
            m_canopy_cv = mean(canopy_cover, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(across(everything(), ~replace(. , is.nan(.), 0)))

# summarizing for grids
grid_sum <- plot_sum %>% 
  filter(pitfall == "Grid") %>% 
  group_by(grid_name, village) %>% 
  summarise_at(vars(4:11), ~mean(., na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(grid_name) %>% 
  mutate(grid_name = factor(grid_name))

# transforming to matrix
grid_sum_mat <- grid_sum %>% 
  column_to_rownames("grid_name") %>%
  ungroup() %>% 
  select(-village) %>%
  as.matrix() %>% 
  decostand(method = "standardize")

#write_csv(grid_sum, "data/data_processed/village_attributes.csv")

## calculating *dis-similarity* between grids
distmat_grid <- as.matrix(vegdist(grid_sum_mat, method = "euclidean"))

# long format
# removing duplicated values
distmat_grid2 <- distmat_grid
distmat_grid2[upper.tri(distmat_grid2)] <- NA
diag(distmat_grid2) <- NA
grid_disimilarity <- melt(distmat_grid2) %>% 
  rename(grid1=Var1, grid2=Var2, grid_attr=value) %>% 
  filter(!(is.na(grid_attr))) %>% 
  arrange(grid1, grid2)

return(grid_disimilarity)
}


##### grid distance #######################################################################

fun_grid_distance <- function(dat) {
  
  dat %<>% arrange(grid)
  grid_dist <- geosphere::distm(dat[3:4], fun = distHaversine)
    rownames(grid_dist) <- dat$grid
    colnames(grid_dist) <- dat$grid
  
  # long format
  # removing duplicated values
  grid_dist2 <- grid_dist
  grid_dist2[upper.tri(grid_dist2)] <- NA
  diag(grid_dist2) <- NA
  grid_distance <- melt(grid_dist2) %>% 
    rename(grid1=Var1, grid2=Var2, grid_dist=value) %>% 
    filter(!(is.na(grid_dist))) %>% 
    arrange(grid1, grid2)
  
  return(grid_distance)
}


##### grid small mammals community #######################################################################

fun_grid_mammals <- function(dat) {
  
  
  # counting the abundance of each host species in every land use and transforming to matrix
  dat_mat <- dat %>% 
    filter(grepl("TMR", animal_id)) %>% 
    rename(host_species = field_identification, grid = habitat_type) %>% 
    filter(host_species != "Rattus rattus") %>% 
    select(host_species, grid) %>% 
    count(grid, host_species) %>% 
    spread(host_species, n, fill = 0) %>% 
    arrange(grid) %>% 
    column_to_rownames("grid") %>% 
    as.matrix()
  
  # calculating *dis-similarity* between grids
  mammals_disimilarity <- as.matrix(vegdist(sqrt(dat_mat), method = "bray"))
    
  # long format
  # removing duplicated values
  mammals_disimilarity2 <- mammals_disimilarity
  mammals_disimilarity2[upper.tri(mammals_disimilarity2)] <- NA
  diag(mammals_disimilarity2) <- NA
  grid_mammals <- melt(mammals_disimilarity2) %>% 
    rename(grid1=Var1, grid2=Var2, sm_community=value) %>% 
    filter(!(is.na(sm_community))) %>% 
    arrange(grid1, grid2)
  
  return(grid_mammals)
}



############################################################################
# main script

# reading the raw data
habitat.raw <- read_csv("data/data_raw/data_small_mammals/Trap_Plots.csv")
plots_location <- read_csv("data/data_raw/data_small_mammals/plots_location.csv")
small_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")

village_names <- unique(plots_location$village)
village_summary <- NULL

# for loop for three villages
for (v in village_names) {
  
  # grid attributes
  grid_attributes <- fun_grid_attributes(habitat.raw %>% filter(village == v))
  # grid distance
  grid_distance <- fun_grid_distance(plots_location %>% filter(village == v))
  # grid community similarity
  grid_mammals <- fun_grid_mammals(small_mammals %>% filter(village == v))
  # combining all variables into one table
  grid_summary <- full_join(grid_attributes, grid_distance, by=c("grid1","grid2")) %>% 
    full_join(grid_mammals, by=c("grid1","grid2")) %>% 
    mutate(village = v) # adding village col
  
  # summary table of three villages together
  village_summary <- rbind(village_summary, grid_summary)
}

# log transformation for the distance variable to match the variables scales
village_summary2 <- village_summary %>% mutate(grid_dist = log(grid_dist))
# saving the results
write_csv(village_summary2, "data/data_processed/village_summary2.csv")


#####
village_summary <- read_csv("data/data_processed/village_summary.csv")
# correlations between variables
psych::pairs.panels(village_summary2 %>% filter(grid1!="village") %>%  select(-village,-grid1,-grid2), ellipses = F, lm = T)


