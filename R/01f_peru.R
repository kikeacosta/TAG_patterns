# source(here("Code", "00_functions.R"))
library(here)
library(readxl)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~
# Peru mortality data
# ~~~~~~~~~~~~~~~~~~~

# Directly from the ministry of health
# https://www.datosabiertos.gob.pe/dataset/informaci%C3%B3n-de-fallecidos-del-sistema-inform%C3%A1tico-nacional-de-defunciones-sinadef-ministerio
# https://www.minsa.gob.pe/reunis/data/defunciones_registradas.asp
# load data
# db_pe <- 
#   read_delim("https://cloud.minsa.gob.pe/s/nqF2irNbFomCLaa/download", delim = "|")
db_pe <- 
  read.csv("Data/peru/SINADEF_DATOS_ABIERTOS.csv")

# data wrangling
db_pe2 <- 
  db_pe %>% 
  select(Sex = SEXO,
         Age = EDAD, 
         Year = AÑO,
         unit_age = TIEMPO.EDAD) %>% 
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
  filter(Year >= 2017 & Year <= 2021) %>% 
  select(Year, Sex, Age, Deaths)

# all sex and all age Peru and Mexico
db_pe_all_ages <- 
  db_pe3 %>% 
  group_by(Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  filter(Sex != "UNK")

db_pe_all_sex <- 
  db_pe3 %>% 
  group_by(Year, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>% 
  filter(Age != "UNK")

db_pe_all_sex_age <- 
  db_pe3 %>% 
  group_by(Year) %>% 
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
  group_by(Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>%
  ungroup()

db_pe5 <- 
  db_pe4 %>% 
  mutate(Country = "Peru", 
         Code = "PER", 
         Source = "peru_sinadef",
         Age = Age %>% as.double()) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Year, Sex, Age)

# saving annual deaths in Peru
write_csv(db_pe5, "Output/peru.csv")

