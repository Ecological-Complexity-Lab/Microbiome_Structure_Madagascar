# script for making a table for ML 

library(tidyverse)
library(magrittr)


rm(list=ls())

# reading the data
# small mammals (hosts)
data_mammals <- read_csv("data/data_raw/data_small_mammals/Terrestrial_Mammals.csv")
data_mammals %<>% mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>% 
  select(host_ID, elevation.obs, sex, mass, age_repro) %>% 
  mutate(sex = as.factor(sex), age_repro = as.factor(age_repro)) 



data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv")
data_asv_filtered <- data_asv %>% 
  filter(host_species == "Rattus rattus" & village == "Mandena" & grid != "village") %>% 
  select(host_ID, asv_ID, grid, season) %>% 
  left_join(data_mammals, by="host_ID") 


    
    
    # ------------
    # making the full table
    
    data_asv_mat <- data_asv_filtered %>% 
      mutate(link = 1) %>% 
      spread(asv_ID, link, fill = 0)  %>% 
      gather("asv_ID","Link", starts_with("ASV"))
  

# saving the final table as .csv
write_csv(data_asv_mat, "data/data_processed/microbiome/ML_microgale_mandena.csv")

