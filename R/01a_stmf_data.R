library(here)
library(tidyverse)
library(countrycode)
source("R/00_functions.R")

# Loading STMF data
# ~~~~~~~~~~~~~~~~~

# downloading the last version of STMF Mortality input data zip 
# this version as of 25 May 2021
download.file("https://www.mortality.org/Public/STMF/Inputs/STMFinput.zip", here("Data/STMFinput.zip"))

# list of country codes in STMF
zipdf <- unzip(here("Data", "STMFinput.zip"), list = TRUE)

# loading all cause deaths from all countries in STMF
db_d <- tibble()
for(i in 1:length(zipdf$Name)){
  csv_file <- zipdf$Name[i]
  print(csv_file)
  temp <- read_csv(unz(here("Data", "STMFinput.zip"), csv_file))
  db_d <- db_d %>% 
    bind_rows(temp)
}

unique(db_d$PopCode)

# countries with full 2020
cts_2020 <- 
  db_d %>% 
  drop_na(Week) %>% 
  filter(Year == 2020) %>% 
  group_by(PopCode) %>% 
  filter(max(Week) >= 52) %>% 
  pull(PopCode) %>% unique()

# countries with full 2021
cts_2021 <- 
  db_d %>% 
  filter(Year == 2021) %>% 
  group_by(PopCode) %>% 
  filter(max(Week) == 52) %>% 
  pull(PopCode) %>% unique()

# filtering periods 2015-2021 with full annual information
dts <- 
  db_d %>% 
  mutate(PopCode = ifelse(PopCode == "a", "NOR", PopCode)) %>% 
  filter(Year %in% 2015:2021) %>% 
  filter(PopCode %in% cts_2020) %>% 
  filter(Year <= 2020 | PopCode %in% cts_2021) %>% 
  select(-Access, -Type, -Area)

unique(dts$PopCode)

# countries with changes in age groups
unique_ages_year <- 
  dts %>% 
  select(PopCode, Age, AgeInterval) %>% 
  unique() %>% 
  group_by(PopCode, Age) %>% 
  summarise(n = n()) %>% 
  filter(n >= 2)

# countries with changes in sex
unique_sex_year <- 
  dts %>% 
  select(PopCode, Sex, Week, Year) %>% 
  unique() %>% 
  group_by(PopCode, Week, Year) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(PopCode, Year, n) %>% 
  unique() %>% 
  filter(n < 3)

rus_sext_19_20 <- 
  dts %>% 
  filter(PopCode == "RUS",
         Sex != "b",
         Year >= 2019) %>% 
  group_by(PopCode, Year, Week, Age, AgeInterval) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "b")

rus <- 
  dts %>% 
  filter(PopCode == "RUS",
         !(Sex == "b" & Year >= 2019)) %>% 
  bind_rows(rus_sext_19_20)

# manual homogenization ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# France: group 0-4 and close at 90
# Russia: group 0-4 and close at 90; calculate both sexes in 2019 and 2020
# Scotland: exclude it
# England and Wales: exclude it
# Northern Ireland: exclude it
# USA: exclude it
exc <- c("USA", "GBR_NIR", "GBRTENW", "GBR_SCO")

dts2 <- 
  dts %>% 
  filter(!PopCode %in% c(exc, "RUS")) %>% 
  bind_rows(rus) %>% 
  mutate(Age = case_when(PopCode %in% c("RUS", "FRATNP") & Age == "1" ~ "0", 
                         PopCode %in% c("RUS", "FRATNP") & Age == "95" ~ "90", 
                         TRUE ~ Age)) %>% 
  group_by(PopCode, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  filter(Age != "UNK") %>% 
  mutate(Sex = recode(Sex,
                      "b" = "t"),
         PopCode = recode(PopCode,
                          "AUS2" = "AUS",
                          "DEUTNP" = "DEU",
                          "FRATNP" = "FRA",
                          "NZL_NP" = "NZL",
                          "GBR_NIR" = "GBR-NIR", 
                          "GBR_SCO" = "GBR-SCO", 
                          "GBRTENW" = "GBR-ENW")) %>%
  rename(Code = PopCode) %>% 
  drop_na(Deaths) %>% 
  mutate(Source = "stmf",
         Country = countrycode(sourcevar = Code, 
                               origin = "iso3c", 
                               destination = "country.name"),
         Country = case_when(Code == "GBR-NIR" ~ "Northern Ireland", 
                             Code == "GBR-SCO" ~ "Scotland", 
                             Code == "GBR-ENW" ~ "England and Wales",
                             Code == "USA" ~ "USA", 
                             TRUE ~ Country)) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Code, Year, Sex, Age)

# saving data
write_csv(dts2, "Output/stmf.csv")


