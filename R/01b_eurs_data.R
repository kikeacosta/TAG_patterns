library(eurostat)
library(tidyverse)
library(countrycode)
IN <- get_eurostat("demo_r_mwk_05")

db_eurs <-
  IN %>% 
  mutate(time = as.character(time)) %>% 
  separate(time, sep = "W", into = c("year","week"), convert = TRUE) %>% 
  dplyr::filter(year >= 2016, 
                year < 2021, 
                age != "UNK") %>% 
  group_by(geo, year) %>% 
  mutate(complete = any(week == 52)) %>% # should I check for 53, or make it conditional?
  ungroup() %>% 
  dplyr::filter(complete) %>% 
  group_by(geo, sex, year, age) %>% 
  summarize(deaths = sum(values), .groups = "drop") %>% 
  mutate(age = recode(age,
                      "TOTAL" = "TOT",
                      "Y_GE90" = "90",
                      "Y_LT5" = "0",
                      "Y5-9" = "5",
                      "Y10-14" = "10",
                      "Y15-19" = "15",
                      "Y20-24" = "20",
                      "Y25-29" = "25",
                      "Y30-34" = "30",
                      "Y35-39" = "35",
                      "Y40-44" = "40",
                      "Y45-49" = "45",
                      "Y50-54" = "50",
                      "Y55-59" = "55",
                      "Y60-64" = "60",
                      "Y65-69" = "65",
                      "Y70-74" = "70",
                      "Y75-79" = "75",
                      "Y80-84" = "80",
                      "Y85-89" = "85"),
         sex = tolower(sex),
         Country = suppressWarnings(countrycode(sourcevar = geo, 
                               origin = "iso2c", 
                               destination = "country.name")),
         Country = case_when(geo == "EL" ~ "Greece",
                             geo == "UK" ~ "United Kingdom",
                             TRUE ~ Country),
         Code = countryname(Country, destination = "iso3c"),
         Source = "eurs") %>% 
  dplyr::select(Country, Code, Year = year, Sex = sex, Age = age, Deaths = deaths, Source)
  
write_csv(db_eurs, "Output/eurs.csv")





  

