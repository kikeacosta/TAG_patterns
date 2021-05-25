library(eurostat)
library(tidyverse)
IN <- get_eurostat("demo_r_mwk_05")

IN$geo %>% unique()

A <-
IN %>% 
  mutate(time = as.character(time)) %>% 
  separate(time, sep = "W", into = c("year","week"), convert = TRUE)
A %>% 
  dplyr::filter(year == 2020) %>% 
  dplyr::pull(age) %>% table()
A %>% 
  dplyr::filter(year == 2020) %>% 
  filter(geo == "GE") %>% 
  dplyr::pull(week) %>% table()

# detect complete years, filter
# then aggregate to years.
B <- 
  A %>% 
  dplyr::filter(year >= 2016, year < 2021) %>% 
  group_by(geo, year) %>% 
  mutate(complete = any(week == 52)) %>% 
  ungroup() %>% 
  dplyr::filter(complete) %>% 
  group_by(geo, sex, year, age) %>% 
  summarize(deaths = sum(values), .groups = "drop")
  
EUR_2020_save <- 
  B %>% 
  dplyr::filter(year == 2020)

D <-
  B %>% 
  dplyr::filter(year != 2020) %>% 
  group_by(geo, sex, age) %>% 
  summarize(deaths = mean(deaths),
            year = mean(year),
            .groups = "drop") %>% 
  bind_rows(EUR_2020_save)






  

