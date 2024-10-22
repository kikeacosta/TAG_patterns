library(readr)
library(tidyverse)
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
         lmx_spinoff = log(dx/Nx))

WM_estimates<- read_excel("Output/EstimatesBySexAge.xlsx", sheet = 2, skip = 5)

WM_pred <- 
  WM_estimates %>% 
  dplyr::filter(measure == "mx",
                `source year` == "Predicted 2020") %>% 
  mutate(lmx_WM = log(mean)) %>% 
  select(iso3, sex, age, lmx_WM)

Data_2019_3 <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_spinoff, lmx_WM) %>% 
  drop_na() %>% 
  gather(lmx_spinoff, lmx_WM, key = source, value = lmx)

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
    geom_line(aes(age, lmx, col = source, linetype = sex), size=1)+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()+
    labs(x = "age", y = "log-mortality")
    
    print(p)
}
dev.off()

test <- 
  Data_2019_3 %>% 
  group_by(iso3, sex, source) %>% 
  summarise(lmx_tot = sum(lmx))


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


diffs <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_spinoff, lmx_WM) %>% 
  drop_na() %>% 
  mutate(diff_lx = exp(lmx_spinoff - lmx_WM))

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
    geom_line(aes(age, diff_lx, linetype = sex),size=1)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_y_log10()+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()+
    labs(x = "age", y = "Incidence rate ratio")
  
  print(p)
}
dev.off()

sex_ratio <- 
  Data_2019_2 %>% 
  left_join(WM_pred) %>% 
  select(Country, iso3, cluster, sex, age, lmx_spinoff, lmx_WM) %>% 
  drop_na() %>% 
  gather(lmx_spinoff, lmx_WM, key = source, value = lmx) %>% 
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
    geom_line(aes(age, sex_ratio, col = source), size=1)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_y_log10()+
    labs(title = paste("Cluster",cl,ctry)) +
    theme_bw()+
    labs(x = "age", y = "Sex ratio")
  
  print(p)
}
dev.off()

