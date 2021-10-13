library(tidyverse)
library(readr)
library(countrycode)
library(readxl)



Deltas1 <- 
  read_delim("Output/Deltas_GC.txt", delim = " ", 
             col_names = c("id", "age","sex","cluster","delta.hat","delta.hatU", "delta.hatL"),
             skip = 1) %>% 
  select(-id)


Data_2019 <- 
  read_delim("Output/data2019_GC.txt", delim = " ", 
             col_names = c("id", "Country", "iso3", "sex", "age", "Nx", "deaths", "cluster", "lmx", "lmx.hat", "lmx.hatU","lmx.hatL"),
             skip = 1) %>% 
  select(-id)

Data_2020 <- 
  read_delim("Output/data2020_GC.txt", delim = " ", 
             col_names = c("id", "Country", "iso3", "sex", "age", "Nx", "deaths", "cluster", "lmx", "lmx.hat", "lmx.hatU","lmx.hatL"),
             skip = 1) %>% 
  select(-id)

WMout <- read_csv("Output/WM_constraints_and_clusters.csv")


WMout2 <- 
  WMout %>% 
  gather(Female, Male, key = sex, value = d_tot) %>% 
  select(iso3, sex, d_tot)

Data_2019_2 <- 
  Data_2019 %>% 
  left_join(Deltas1) %>% 
  mutate(lmx_cluster = lmx + delta.hat,
         lmx_clusterU = lmx + delta.hatU,
         lmx_clusterL = lmx + delta.hatL,
         dx_cluster = exp(lmx_cluster) * Nx,
         dx_clusterU = exp(lmx_clusterU) * Nx,
         dx_clusterL = exp(lmx_clusterL) * Nx) %>% 
  group_by(Country, iso3, cluster, sex) %>% 
  mutate(px = dx_cluster / sum(dx_cluster)) %>% 
  ungroup() %>% 
  left_join(WMout2) %>% 
  mutate(dx = px * d_tot,
         lmx_scaled = log(dx/Nx))

WM_estimates<- read_excel("Data/EstimatesBySexAge.xlsx", sheet = 2, skip = 5)

WM_pred <- 
  WM_estimates %>% 
  dplyr::filter(measure == "mx",
                `source year` == "Predicted 2020") %>% 
  mutate(lmx_pred_w = log(mean)) %>% 
  select(iso3, sex, age, lmx_pred_w)

Data_2019_3 <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_scaled, lmx_pred_w) %>% 
  drop_na() %>% 
  gather(lmx_scaled, lmx_pred_w, key = source, value = lmx)

write_csv(Data_2019_3, "Output/age_sex_compare3.csv")

pdf("Figures/in_sample_compare2.pdf")

isos <- unique(Data_2019_3$iso3)
for(is in isos){
  
  chunk <- 
    Data_2019_3 %>% 
    filter(iso3 == is)
  
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl  <-  chunk$cluster[1]
  p <- 
    chunk %>% 
    ggplot()+
    geom_line(aes(age, lmx, col = source, linetype = sex))+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()
  
  print(p)
}
dev.off()

test <- 
  Data_2019_3 %>% 
  group_by(iso3, sex, source) %>% 
  summarise(lmx_tot = sum(lmx))



# verifying the scaling is working well with William's estimates
WM_deaths <- 
  WM_estimates %>% 
  dplyr::filter(measure == "deaths",
                `source year` == "Predicted 2020") %>% 
  select(iso3, sex, age, deaths = mean) %>% 
  group_by(iso3, sex) %>% 
  summarise(d_tot_wm = sum(deaths))


own_deaths <- 
  Data_2019_2 %>% 
  select(iso3, sex, age, dx) %>% 
  group_by(iso3, sex) %>% 
  summarise(d_tot_gc = sum(dx))

deaths <- 
  own_deaths %>% 
  left_join(WM_deaths) %>% 
  drop_na() %>% 
  mutate(diff = d_tot_gc - d_tot_wm)
# All good!!

# differnces of estimates
diffs <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_scaled, lmx_pred_w) %>% 
  drop_na() %>% 
  mutate(diff_lx = exp(lmx_scaled - lmx_pred_w))

pdf("Figures/all_diffs.pdf")

isos <- unique(diffs$iso3)
for(is in isos){
  
  chunk <- 
    diffs %>% 
    filter(iso3 == is)
  
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl  <-  chunk$cluster[1]
  p <- 
    chunk %>% 
    ggplot()+
    geom_line(aes(age, diff_lx, linetype = sex))+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_y_log10()+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()
  
  print(p)
}
dev.off()

# sex ratios
sex_ratio <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_scaled, lmx_pred_w) %>% 
  drop_na() %>% 
  gather(lmx_scaled, lmx_pred_w, key = source, value = lmx) %>% 
  spread(sex, lmx) %>% 
  mutate(sex_ratio = exp(Male - Female))

pdf("Figures/all_sex_ratios.pdf")

isos <- unique(sex_ratio$iso3)
for(is in isos){
  
  chunk <- 
    sex_ratio %>% 
    filter(iso3 == is)
  
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl  <-  chunk$cluster[1]
  p <- 
    chunk %>% 
    ggplot()+
    geom_line(aes(age, sex_ratio, linetype = source))+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_y_log10()+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()
  
  print(p)
}
dev.off()

# =============================
# saving estimates in dx and Nx
WM_dx <- 
  WM_estimates %>% 
  dplyr::filter(measure == "deaths",
                `source year` %in% c("Predicted 2020", "Expected 2020", "Observed 2020")) %>% 
  select(Country, iso3, sex, age, Nx, measure = 4, dx = mean) %>% 
  mutate(measure = recode(measure,
                          "Predicted 2020" = "predicted_dx_tag1",
                          "Expected 2020" = "baseline_dx_tag1", 
                          "Observed 2020" = "observed")) %>% 
  spread(measure, dx)

dx_compare <- 
  Data_2019_2 %>% 
  select(Country, iso3, sex, age, predicted_dx_tag2 = dx) %>% 
  left_join(WM_dx) %>% 
  drop_na(baseline_dx_tag1) %>% 
  drop_na(predicted_dx_tag1)

write_csv(dx_compare, "Output/age_sex_deaths_comparison.csv")


