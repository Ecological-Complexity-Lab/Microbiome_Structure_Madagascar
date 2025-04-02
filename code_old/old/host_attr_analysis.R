# host attributes analysis
# exploring the attributes of hosts characterized by different core/non-core microbes

library(tidyverse)
library(magrittr)
library(ggplot2)

rm(list=ls())


# reading microbiome data
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv") %>% 
  filter(host_species == "Rattus rattus" & grid!="village" & village == "Mandena")

# calculating ASV degree for each village
asv_degree <- data_asv %>% 
  count( asv_ID) %>% rename(asv_degree=n)
data_asv %<>% left_join(asv_degree, by="asv_ID")
# modules
data_modules <- fun_modularity_analysis(data_asv)

n_landuse <- data_modules %>% 
  group_by(host_group) %>% 
  summarise(n_grid = n_distinct(grid))

host_grid <- data_modules %>% distinct(host_ID, host_group)

# calculating ASV degree for each village
asv_degree <- data_asv %>% 
  count(village, asv_ID)  

# calculating the median ASV degree for each host
host_median_degree <- data_asv %>% 
  group_by(host_ID) %>% 
  summarise(median_degree = median(asv_degree), mean_degree = mean(asv_degree))

# adding host attributes
host_attr <- read_csv("data/data_processed/small_mammals/small_mammals_attributes.csv")

host_final_data <- host_median_degree %>% 
  left_join(host_attr, by="host_ID") %>% 
  left_join(host_grid, by="host_ID") %>% 
  left_join(n_landuse, by="host_group") %>% 
  #mutate(age_repro = factor(age_repro)) %>% 
  mutate(pre = ifelse(n_grid<3,"low", ifelse(n_grid>4,"high", "medium"))) %>% 
  mutate(pre = factor(pre, levels = c("low","medium","high"))) %>%   
  mutate(sex2 = ifelse(sex=="male",1,0)) %>% 
  mutate(age = ifelse(age_repro==1,1,0))

# exploring the distribution
a=quantile(host_final_data$median_degree, c(0.1,0.5,0.9))

# difference between groups
diff_results <- host_final_data %>%
  #filter(sex == "male") %>% 
  group_by(pre) %>% 
  summarise(n = mean(age_repro))
count(age)

host_final_data %>% 
ggplot(aes(x=pre, y=sex2)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'), title = element_text(size = 20), legend.position = "none")

# logistic regression
library(lme4)
sex_reg <- glmer(median_degree ~ age_repro + sex + (1|village) + (1|grid), data = host_final_data, family = "gaussian")
summary(sex_reg)

library(nlme)
sex_reg <- lme(median_degree ~ age + sex, random=~1|village , data = host_final_data, method = "ML", na.action = na.fail)
summary(sex_reg)
