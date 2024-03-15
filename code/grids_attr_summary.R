# script for analyze grids attributes
# *taken from Kayla script*

library(tidyverse)
library(magrittr)
library(vegan)

setwd("GitHub/Small_Mammals_Microbiome")
rm(list=ls())

############################################################################


habitat.raw <- read_csv("data_raw/Trap_Plots.csv")

habitat <- habitat.raw %>% 
  filter(village == "Mandena") %>% 
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
  select(village, grid_name, pitfall, plot, square, tree_height, tree_dbh)

## create a new df individual tree measurements
squares <- habitat %>% 
  mutate(tree = ifelse(is.na(tree_dbh), 0, 1)) %>% ## count the number of trees in each square
  group_by(village, grid_name, plot, square) %>% 
  mutate(n_trees = sum(tree)) %>% 
  ungroup() %>% 
  select(-tree_height, -tree_dbh, -tree) %>% # remove tree info from habitat
  distinct()

plot_sum <- squares %>% 
  mutate(across(c(n_trees, dead_logs, 
                  liana_number, herbaceous_height, 
                  herbaceous_cover, canopy_cover))) %>% 
  group_by(pitfall, pitfall_line, village, grid_name, plot) %>% 
  summarise(n_trees = sum(n_trees, na.rm=T),
            n_logs = sum(dead_logs, na.rm=T),
            n_liana = sum(liana_number, na.rm=T),
            m_herb_ht = mean(herbaceous_height, na.rm=T),
            m_herb_cv = mean(herbaceous_cover, na.rm = T),
            m_canopy_cv = mean(canopy_cover, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(across(everything(), ~replace(. , is.nan(.), 0)))


grid_sum <- plot_sum %>% 
  filter(pitfall == "Grid") %>% 
  group_by(grid_name, village) %>% 
  summarise_at(vars(4:9), ~mean(., na.rm=TRUE)) %>% 
  ungroup() %>% 
  #add_row(grid_name="village",village="Mandena",n_trees=0,n_logs=0,n_liana=0,m_herb_ht=0,m_herb_cv=0,m_canopy_cv=0) %>% 
  arrange(grid_name) %>% 
  mutate(grid_name = factor(grid_name))
#write_csv(grid_sum, "data_processed/grid_attr_summary.csv")

grid_sum_mat <- grid_sum %>% 
  #filter(pitfall == "Grid") %>% 
  #column_to_rownames("grid_name") %>%
  ungroup() %>% 
  select(-grid_name,-village) %>%
  as.matrix()

pitfall_sum_mat <- grid_sum %>% 
  filter(pitfall == "Pitfall") %>% 
  column_to_rownames("grid_name") %>%
  select(-pitfall) %>%
  as.matrix()

## dist matrix
distmat_grid <- as.matrix(vegdist(grid_sum_mat, method = "bray")) 
distmat_pitfall <- 1-vegdist(sqrt(pitfall_sum_mat), method = "bray") 


## nmds
nms <- metaMDS(distmat_grid,
               k = 3, 
               maxit = 999,
               trymax = 500,
               wascores = T)
## plot ordination grids
ordiplot(nms, type = "none")
points(nms, display = "sites", bg = grid_sum$grid_name)
text(nms, display = "sites", bg = grid_sum$grid_name)




# mantel test
mantel_grid_attr <- ecodist::mantel(distmat_grid ~ as.dist(grid_dist[-7,-7]))
mantel_grid_attr
mantel_grid_attr <- ecodist::mantel(as.dist(grid_modules_dist[-7,-7]) ~ distmat_grid)
mantel_grid_attr
mantel_grid_attr <- ecodist::mantel(as.dist(grid_modules_dist) ~ as.dist(grid_dist))
mantel_grid_attr
mmm <- ecodist::MRM(as.dist(grid_modules_dist[-7,-7]) ~ as.dist(grid_dist[-7,-7]) + distmat_grid)
mmm
