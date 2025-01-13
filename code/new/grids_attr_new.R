# script for analyze grids attributes
# *taken from Kayla script*

library(tidyverse)
library(magrittr)
library(vegan)
library(geosphere)
library(reshape2)
library(psych)
library(ggrepel)

rm(list=ls())


##### grid vegetation #######################################################################

habitat.raw <- read_csv("data/data_raw/data_small_mammals/Trap_Plots.csv")

habitat <- habitat.raw %>% 
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
  select(village, grid_name, season, plot, square, dead_logs, tree_dbh, tree_height, #omiting crown distance
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
  select(village, grid_name, season, pitfall, plot, square, tree_height, tree_dbh) %>% 
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

# name index
grid_village_name <- grid_sum %>% 
  select(village, grid_name) %>% 
  unite(col="grid_village", village:grid_name, remove = FALSE)

# transforming to matrix
grid_sum_mat <- grid_sum %>% 
  unite(col="grid_village", c("village","grid_name"), remove = TRUE) %>%
  column_to_rownames("grid_village") %>%
  ungroup() %>% 
  as.matrix() 

# PCA
# variables transformed to centered 0 and rescale to have unit variance before the analysis 
pca_veg <- prcomp(grid_sum_mat, center = TRUE, scale. = TRUE)

vegetation_pca <- as.data.frame(pca_veg$x) %>% 
  rownames_to_column("grid_village") %>% 
  left_join(grid_village_name, by="grid_village") %>% 
  select(grid_village,village,grid_name,PC1,PC2) %>% 
  dplyr::rename(grid=grid_name) %>% 
  rename_with(~paste("veg",.,sep="_"), .cols=starts_with("PC"))

# explained variance
ex_var <- pca_veg$sdev ^2 
prop_ex_var <- ex_var/sum(ex_var)*100
prop_ex_var

pca_data <- data.frame(pca_veg$x[, 1:2],  # Select the first two principal components
                       village = as.factor(grid_sum$village), 
                       land_use = as.factor(grid_sum$grid_name)) 

loadings <- data.frame(pca_veg$rotation[, 1:2], variable = rownames(pca_veg$rotation))

#pdf(file = 'results/figure_S3.pdf')
veg_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 5, alpha = 0.7, aes(fill = land_use, shape = village)) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*6, yend = PC2*6), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +  # Vectors as arrows
  geom_text(data = loadings, aes(x = PC1 * 7, y = PC2 * 7, label = variable),
            color = "black", size = 4) +
  labs(x = paste("PC1 (",round(prop_ex_var[1],2),"%)", sep = ""),
       y = paste("PC2 (",round(prop_ex_var[2],2),"%)", sep = ""),
       fill = "Land-use", shape="Village") +
  scale_fill_manual(values = c("red", "blue", "green", "orange", "purple","yellow", "black", "grey")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black"), 
        panel.grid = element_blank(), aspect.ratio = 1, title = element_text(size = 14))
#dev.off()
veg_plot


#####
# small mammals

small_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")

# counting the abundance of each host species in every land use and transforming to matrix
dat_mat <- small_mammals %>% 
  filter(grepl("TMR", animal_id)) %>% 
  dplyr::rename(host_species = field_identification, grid = habitat_type) %>% 
  filter(host_species != "Rattus rattus") %>% 
  select(host_species,village, grid) %>% 
  count(village, grid, host_species) %>% 
  spread(host_species, n, fill = 0) 

# name index
grid_village_name <- dat_mat %>% 
  select(village, grid) %>% 
  unite(col="grid_village", village:grid, remove = FALSE)

# final matrix
grid_dat_mat <- dat_mat %>% 
  unite(col="grid_village", village:grid, remove = TRUE) %>% 
  column_to_rownames("grid_village") %>% 
  as.matrix()

# PCA
# variables transformed to centered 0 and rescale to have unit variance before the analysis 
pca_sm <- prcomp(sqrt(grid_dat_mat), center = FALSE, scale. = FALSE)

sm_community_pca <- as.data.frame(pca_sm$x) %>% 
  rownames_to_column("grid_village") %>% 
  left_join(grid_village_name, by="grid_village") %>% 
  select(grid_village,PC1,PC2) %>% 
  rename_with(~paste("sm",.,sep="_"), .cols=starts_with("PC"))

# explained variance
ex_var <- pca_sm$sdev ^2 
prop_ex_var <- ex_var/sum(ex_var)*100
prop_ex_var



##### grid distance from village #######################################################################

fun_distance_village <- function(dat) {
  # calculating the difference between grids in the distance to village
  
  dat %<>% arrange(grid)
  vil = dat %>% filter(grid=="village") %>% select(longitude, latitude)
  dist_village_grids <- as.data.frame(geosphere::distm(dat[3:4], vil , fun = distHaversine))
  rownames(dist_village_grids) <- dat$grid
  colnames(dist_village_grids) <- "dist_to_village"
  dist_village_grids <- rownames_to_column(dist_village_grids, "grid") %>% 
    mutate(village = dat$village[1])
  
  return(dist_village_grids)
}

plots_location <- read_csv("data/data_raw/data_small_mammals/plots_location.csv")
village_names <- sort(unique(plots_location$village))
village_dist <- NULL
# for loop for three villages
for (v in village_names) {
  
  # distance from village
  village_distance <- fun_distance_village(plots_location %>% filter(village == v))
  village_dist <- rbind(village_dist, village_distance)
}

village_distance <- village_dist %>% 
  unite(col="grid_village", c(village,grid), remove = TRUE)


##### distance between grids #######################################################################

# calculating pairwise distance between all grids from all villages
dist_grids <- as.data.frame(geosphere::distm(plots_location[3:4], fun = distHaversine))
site_names <- plots_location %>% unite(col="site", c(village,grid), remove = TRUE)
rownames(dist_grids) <-  site_names$site
mds_distance <- cmdscale(dist_grids, k = 1) %>% 
  as.data.frame() %>% 
  rename(distance = V1) %>% 
  rownames_to_column("grid_village")



##### grid elevation #######################################################################
  
  # finding elevation for each grid (average value of captured animals)
  grid_elevation <- small_mammals %>% 
    filter(grepl("TMR", animal_id)) %>% 
    dplyr::rename(host_species = field_identification, grid = habitat_type) %>% 
    group_by(village, grid) %>% 
    summarise(elevation = mean(elevation.obs)) %>% 
  unite(col="grid_village", c(village,grid), remove = TRUE)



############################################################################
# combine variables

village_summary <- vegetation_pca %>% 
  full_join(sm_community_pca, by="grid_village") %>% 
  full_join(village_distance, by="grid_village") %>% 
  full_join(grid_elevation, by="grid_village") %>% 
  full_join(mds_distance, by="grid_village")

# saving the results
write_csv(village_summary, "data/data_processed/village_summary_new.csv")



#####
village_summary <- read_csv("data/data_processed/village_summary.csv")
# correlations between variables
psych::pairs.panels(village_summary %>%  select(-grid_village), ellipses = F, lm = F)


