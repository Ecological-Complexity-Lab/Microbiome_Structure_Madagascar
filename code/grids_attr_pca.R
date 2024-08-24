# script for calculating grid attributes

rm(list=ls())

#####
# vegetation

habitat.raw <- read_csv("data/data_raw/data_small_mammals/Trap_Plots.csv") %>% 
  mutate(season = ifelse(is.na(season), 3, season))

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
  group_by(village, grid_name,season, plot, square) %>%
  summarise(tree_height = mean(tree_height, na.rm=T),
            tree_dbh = mean(tree_dbh, na.rm=T))

## create a new df individual tree measurements
squares <- habitat %>% 
  mutate(tree = ifelse(is.na(tree_dbh), 0, 1)) %>% ## count the number of trees in each square
  group_by(village, grid_name, season, plot, square) %>% 
  mutate(n_trees = sum(tree)) %>% 
  ungroup() %>% 
  select(-tree_height, -tree_dbh, -tree) %>% # remove tree info from habitat
  distinct() %>% 
  left_join(trees, by=c("village","grid_name","season","plot","square"))

# summarizing for plots
plot_sum <- squares %>% 
  mutate(across(c(n_trees, tree_height, tree_dbh, dead_logs, 
                  liana_number, herbaceous_height, 
                  herbaceous_cover, canopy_cover))) %>% 
  group_by(pitfall, pitfall_line, village, grid_name, season, plot) %>% 
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
  group_by(grid_name, season, village) %>% 
  summarise_at(vars(4:11), ~mean(., na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(grid_name) %>% 
  mutate(grid_name = factor(grid_name))

# name index
grid_season_name <- grid_sum %>% 
  select(grid_name, season, village) %>% 
  unite(col="grid_season", grid_name:village, remove = FALSE)

# transforming to matrix
grid_sum_mat <- grid_sum %>% 
  unite(col="grid_season", grid_name:village, remove = TRUE) %>%
  column_to_rownames("grid_season") %>%
  ungroup() %>% 
  as.matrix() 

# PCA
# variables transformed to centered 0 and rescale to have unit variance before the analysis 
pca_veg <- prcomp(grid_sum_mat, center = TRUE, scale. = TRUE)

vegetation_pca <- as.data.frame(pca_veg$x) %>% 
  rownames_to_column("grid_season") %>% 
  left_join(grid_season_name, by="grid_season")

# explained variance
ex_var <- pca_veg$sdev ^2 
prop_ex_var <- ex_var/sum(ex_var)*100
prop_ex_var

write_csv(vegetation_pca, "data/data_processed/pca_vegetation_3_villages_by_grid_season.csv")

# plot
library(ggplot2)

pca_data <- data.frame(pca_veg$x[, 1:2],  # Select the first two principal components
                       village = as.factor(grid_sum$village),  # Replace with your actual village factor
                       season = as.factor(grid_sum$season),    # Replace with your actual season factor
                       land_use = as.factor(grid_sum$grid_name)) # Replace with your actual land use factor

loadings <- data.frame(pca_veg$rotation[, 1:2], variable = rownames(pca_veg$rotation))

scaling_factor <- 7

veg_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, aes(color = land_use, shape = village, fill = season)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*scaling_factor, yend = PC2*scaling_factor), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +  # Vectors as arrows
  geom_text(data = loadings, aes(x = PC1 * scaling_factor* 1.2, y = PC2 * scaling_factor * 1.2, label = variable),
            color = "black", size = 4) +  # Variable names near the arrows
  labs(title = "PCA vegetation (per land use, village and season)",
       x = paste("PC1 (",round(prop_ex_var[1],2),"% variance explained)", sep = ""),
       y = paste("PC2 (",round(prop_ex_var[2],2),"% variance explained)", sep = "")) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "black", "grey")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("white", "black", "grey"))  +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  stat_ellipse(geom = "polygon",
               aes(color = land_use), # use inherit.aes=F if needed
               alpha=0.1) + # Add ellipses by land_use
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
saveRDS(veg_plot, "pca_veg_plot.Rds")


#####
# small mammals

small_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")

# counting the abundance of each host species in every land use and transforming to matrix
dat_mat <- small_mammals %>% 
  filter(grepl("TMR", animal_id)) %>% 
  rename(host_species = field_identification, grid = habitat_type) %>% 
  filter(host_species != "Rattus rattus") %>% 
  select(host_species,village, grid, season) %>% 
  count(village, grid, season, host_species) %>% 
  spread(host_species, n, fill = 0) 
  
# name index
grid_season_name <- dat_mat %>% 
  select(village, grid, season) %>% 
  unite(col="grid_name", village:season, remove = FALSE)

# final matrix
grid_dat_mat <- dat_mat %>% 
  unite(col="grid_name", village:season, remove = TRUE) %>% 
  column_to_rownames("grid_name") %>% 
  as.matrix()

# PCA
# variables transformed to centered 0 and rescale to have unit variance before the analysis 
pca_sm <- prcomp(grid_dat_mat, center = TRUE, scale. = TRUE)

sm_community_pca <- as.data.frame(pca_sm$x) %>% 
  rownames_to_column("grid_name") %>% 
  left_join(grid_season_name, by="grid_name")

# explained variance
ex_var <- pca_sm$sdev ^2 
prop_ex_var <- ex_var/sum(ex_var)*100
prop_ex_var

write_csv(sm_community_pca, "data/data_processed/pca_sm_community_3_villages_by_grid_season.csv")

# plot
pca_data <- data.frame(pca_sm$x[, 1:2],  # Select the first two principal components
                       village = as.factor(dat_mat$village),  # Replace with your actual village factor
                       season = as.factor(dat_mat$season),    # Replace with your actual season factor
                       land_use = as.factor(dat_mat$grid)) # Replace with your actual land use factor

loadings <- data.frame(pca_sm$rotation[, 1:2], variable = rownames(pca_sm$rotation))

scaling_factor <- 10

sm_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, aes(color = land_use, shape = village, fill = season)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*scaling_factor, yend = PC2*scaling_factor), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +  # Vectors as arrows
  geom_text(data = loadings, aes(x = PC1 * scaling_factor* 1.2, y = PC2 * scaling_factor * 1.2, label = variable),
            color = "black", size = 3) +  # Variable names near the arrows
  labs(title = "PCA vegetation (per land use, village and season)",
       x = paste("PC1 (",round(prop_ex_var[1],2),"% variance explained)", sep = ""),
       y = paste("PC2 (",round(prop_ex_var[2],2),"% variance explained)", sep = "")) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "black", "cyan")) +  # Customize colors as needed
  scale_shape_manual(values = c(21, 22, 23)) +              # Customize shapes as needed
  scale_fill_manual(values = c("white", "black", "grey"))  +  # Customize fills as needed
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  stat_ellipse(geom = "polygon",
               aes(color = land_use),  # use inherit.aes=F if needed
               alpha=0.1) + # Add ellipses by land_use
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

saveRDS(sm_plot, "pca_sm_plot.Rds")


#####
# microbiome

data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.001_p0.01_th5000.csv")

# microbes taxonomy
tax <- read_delim("data/data_raw/data_microbiome/ASVs_taxonomy_new.tsv") %>% 
  dplyr::rename(asv_ID = ASV)

data_asv_tax <- data_asv %>% 
  left_join(tax, by="asv_ID")

# how much classified families?
families <- data_asv_tax %>% 
  distinct(asv_ID, Family) %>% 
  mutate(classified = !is.na(Family))
table(families$classified)  # FALSE= unknown Family
# 2518/2808 = 0.896

# how many families?
length(unique(data_asv_tax$Family))

# aggregating by Family level
data_asv_family <- data_asv_tax %>% 
  group_by(host_ID, Family) %>% 
  summarise(reads = sum(reads)) %>% 
  filter(!is.na(Family))

# matrix
data_asv_family_mat <- data_asv_family %>% 
  spread(Family, reads, fill = 0) %>% 
  column_to_rownames("host_ID") %>% 
  as.matrix()

# PCA
pca_microbiome <- prcomp(data_asv_family_mat, center = FALSE, scale. = FALSE)

microbiome_community_pca <- as.data.frame(pca_microbiome$x) %>% 
  rownames_to_column("host_ID") %>% 
  mutate(host_ID = as.double(host_ID)) %>% 
  left_join(data_asv %>% distinct(host_ID,village,grid,season), by="host_ID")

# explained variance
ex_var <- pca_microbiome$sdev ^2 
prop_ex_var <- ex_var/sum(ex_var)*100
prop_ex_var

loadings <- as.data.frame(pca_microbiome$rotation[, 1:2]) %>% rownames_to_column("Family") %>% 
  mutate(PC1_abs = abs(PC1), PC2_abs = abs(PC2))

loading_top <- loadings %>% 
  slice_max(n = 5, order_by = PC1_abs) %>% 
  bind_rows(loadings %>% slice_max(n = 5, order_by = PC2_abs))

write_csv(sm_community_pca, "data/data_processed/pca_microbiome_community_pca.csv")

# plotting
scaling_factor <- 1.5

ggplot(microbiome_community_pca, aes(x = PC1, y = PC2, color = village)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_segment(data = loading_top, aes(x = 0, y = 0, xend = PC1*scaling_factor, yend = PC2*scaling_factor), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +  
  geom_text(data = loadings, aes(x = PC1 * scaling_factor, y = PC2 * scaling_factor, label = Family),
            color = "black", size = 3) +  
  labs(title = "PCA host microbiome",
       x = paste("PC1 (",round(prop_ex_var[1],2),"% variance explained)", sep = ""),
       y = paste("PC2 (",round(prop_ex_var[2],2),"% variance explained)", sep = "")) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
