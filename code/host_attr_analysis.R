# host attributes analysis
# exploring the attributes of hosts characterized by different core/non-core microbes

library(tidyverse)
library(magrittr)
library(ggplot2)

rm(list=ls())


# reading microbiome data
data_asv <- read_csv("data/data_processed/microbiome/data_asv_rra0.01_th1000.csv") %>% 
  filter(host_species == "Rattus rattus" & grid!="village")

# calculating ASV degree for each village
asv_degree <- data_asv %>% 
  count(village, asv_ID)  

# calculating the median ASV degree for each host
host_median_degree <- data_asv %>% 
  left_join(asv_degree, by=c("village","asv_ID")) %>% 
  group_by(host_ID) %>% 
  summarise(median_degree = median(n), mean_degree = mean(n))

# adding host attributes
host_attr <- read_csv("data/data_processed/small_mammals/small_mammals_attributes.csv")

host_final_data <- host_median_degree %>% 
  left_join(host_attr, by="host_ID") %>% 
  mutate(age_repro = factor(age_repro)) %>% 
  mutate(pre = ifelse(median_degree<5,"low", ifelse(median_degree>27,"high", "medium"))) %>% 
  mutate(pre = factor(pre, levels = c("low","medium","high"))) %>% 
  mutate(sex2 = ifelse(sex=="male",1,0)) %>% 
  mutate(age = ifelse(age_repro==1,1,0))

# exploring the distribution
a=quantile(host_final_data$median_degree, c(0.1,0.5,0.9))

# difference between groups
host_final_data %>% 
ggplot(aes(x=age_repro, y=median_degree, fill=age_repro)) + 
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
