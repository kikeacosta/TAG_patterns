# source(here("Code", "00_functions.R"))
library(here)
library(readxl)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~
# Peru mortality data
# ~~~~~~~~~~~~~~~~~~~

# Directly from the ministry of health
# https://www.datosabiertos.gob.pe/dataset/informaci%C3%B3n-de-fallecidos-del-sistema-inform%C3%A1tico-nacional-de-defunciones-sinadef-ministerio
# load data
db_pe <- 
  read_delim("https://cloud.minsa.gob.pe/s/nqF2irNbFomCLaa/download", delim = "|")

# db_pe <- 
#   read_delim(here("Data", "Peru", "fallecidos_sinadef.csv"), delim = "|")

# data wrangling
db_pe2 <- db_pe %>% 
  select(Sex = SEXO,
         Age = EDAD, 
         Year = AÑO,
         unit_age = 'TIEMPO EDAD') %>% 
  mutate(Sex = recode(Sex,
                      "MASCULINO" = "m",
                      "FEMENINO" = "f"))

db_pe3 <- 
  db_pe2 %>%
  mutate(Sex = case_when(Sex == "m" | Sex == "f" ~ Sex,
                         TRUE ~ "UNK"),
         Age = ifelse(unit_age == "AÑOS", Age, "0"),
         Age = as.integer(Age),
         Age = case_when(Age > 100 ~ "100",
                         is.na(Age) | Age == 999 ~ "UNK", 
                         TRUE ~ as.character(Age)),
         Year = as.double(Year)) %>% 
  group_by(Year, Sex, Age) %>% 
  summarise(Deaths = n()) %>% 
  ungroup() %>% 
  filter(Year >= 2017 & Year <= 2020) %>% 
  mutate(Country = "Peru") %>% 
  select(Country, Year, Sex, Age, Deaths)

# all sex and all age Peru and Mexico
db_pe_all_ages <- 
  db_pe3 %>% 
  group_by(Country, Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  filter(Sex != "UNK")

db_pe_all_sex <- 
  db_pe3 %>% 
  group_by(Country, Year, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>% 
  filter(Age != "UNK")

db_pe_all_sex_age <- 
  db_pe3 %>% 
  group_by(Country, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t",
         Age = "TOT")

db_pe4 <- 
  db_pe3 %>% 
  filter(Age != "UNK",
         Sex != "UNK") %>% 
  bind_rows(db_pe_all_ages,
            db_pe_all_sex,
            db_pe_all_sex_age) %>% 
  arrange(Country, Year, Sex, suppressWarnings(as.numeric(Age))) %>% 
  mutate(Code = "PER", Source = "country_public")

# saving annual deaths in Mexico and Peru
write_csv(db_pe4, "Output/peru.csv")

# db_pe_mx_adj <- 
#   db_pe_mx2 %>% 
#   group_by(Country, Sex, Year) %>% 
#   do(rescale_age(chunk = .data)) %>% 
#   ungroup() %>% 
#   group_by(Country, Age, Year) %>% 
#   do(rescale_sex(chunk = .data)) %>% 
#   ungroup() %>% 
#   arrange(Country, Year, Sex, suppressWarnings(as.numeric(Age)))
# 
# test <- 
#   db_pe_mx_adj %>% 
#   group_by(Country, Year) %>% 
#   summarise(Deaths = sum(Deaths)) %>% 
#   ungroup() %>% 
#   left_join(db_pe_mx %>% 
#               group_by(Country, Year) %>% 
#               summarise(Deaths = sum(Deaths)) %>% 
#               ungroup() %>% 
#               rename(Deaths_org = Deaths)) 
#   
# write_csv(db_pe_mx_adj, "Output/pe_mx_annual_deaths.csv")
# 
