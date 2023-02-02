library(eurostat)
source("R/00_functions.R")
# library(tidyverse)
# library(countrycode)
IN <- get_eurostat("demo_r_mwk_05")

db_eurs <-
  IN %>% 
  mutate(time = as.character(time)) %>% 
  separate(time, sep = "W", into = c("year","week"), convert = TRUE) %>% 
  dplyr::filter(year >= 2015, 
                year <= 2021, 
                age != "UNK") %>% 
  group_by(geo, year) %>% 
  mutate(complete = any(week == 52)) %>% # should I check for 53, or make it conditional?
  ungroup() %>% 
  dplyr::filter(complete) %>% 
  group_by(geo) %>% 
  filter(max(year) >= 2020) %>% 
  ungroup() %>% 
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
  arrange(Country, year, sex, suppressWarnings(as.integer(age))) %>% 
  dplyr::select(Country, Code, Year = year, Sex = sex, Age = age, Deaths = deaths, Source)
  
# re-scaling ages and sexes
db_eurs2 <- 
  db_eurs %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Country, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Code, Year, Sex, Age)

write_csv(db_eurs2, "data_inter/eurs.csv")





  

