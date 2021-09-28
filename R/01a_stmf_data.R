library(here)
library(tidyverse)
library(countrycode)

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
# a PopCode "a"
test <- 
  db_d %>% 
  mutate(id = 1:n())

# it is Norway in week 18, age 90, sex "b"

# temporal fix: Russia has no total sex for 2019 and 2020
test_rus <- 
  db_d %>% 
  filter(PopCode == "RUS",
         Sex == "b",
         Age == "TOT",
         Year >= 2015)

rus_19_20_sex_b <- 
  db_d %>% 
  filter(PopCode == "RUS",
         Sex != "b",
         Year >= 2019) %>% 
  select(-Access, -Type, -AgeInterval, -Area) %>% 
  group_by(PopCode, Year, Week, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "b")

rus <- 
  db_d %>% 
  select(-Access, -Type, -AgeInterval, -Area) %>% 
  filter(PopCode == "RUS",
         Year >= 2015,
         !(Year >= 2019 & Sex == "b")) %>% 
  bind_rows(rus_19_20_sex_b) %>% 
  arrange(Year, Week, Sex)
  
db_d2 <- 
  db_d %>% 
  filter(PopCode != "RUS") %>% 
  select(-Access, -Type, -AgeInterval, -Area) %>% 
  mutate(PopCode = ifelse(PopCode == "a", "NOR", PopCode)) %>% 
  bind_rows(rus)
unique(db_d2$PopCode)

db_stmf <- 
  db_d2 %>% 
  mutate(PopCode = recode(PopCode,
                          "AUS2" = "AUS",
                          "DEUTNP" = "DEU",
                          "FRATNP" = "FRA",
                          "NZL_NP" = "NZL"),
         Sex = recode(Sex,
                      "b" = "t"), 
         Week = as.character(Week)) %>%
  rename(Code = PopCode) %>% 
  drop_na(Deaths) %>% 
  group_by(Code, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  filter(Age != "UNK",
         Year >= 2015,
         Year <= 2020) %>% 
  mutate(Source = "stmf",
         Country = countrycode(sourcevar = Code, 
                               origin = "iso3c", 
                               destination = "country.name"),
         Country = case_when(Code == "GBR_NIR" ~ "Northern Ireland", 
                             Code == "GBR_SCO" ~ "Scotland", 
                             Code == "GBRTENW" ~ "England and Wales",
                             Code == "USA" ~ "USA", 
                             TRUE ~ Country)) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)


write_csv(db_stmf, "Output/stmf.csv")
