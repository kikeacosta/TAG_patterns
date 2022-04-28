library(here)
library(readxl)
source(here("Code", "00_functions.R"))

db <- 
  read_csv("Data/icd_raw.csv",
           col_types = cols(.default = "c"))

cts20 <- 
  db %>% 
  filter(Year == 2020 & Frmat != "09") %>% 
  pull(Country) %>% 
  unique


db20 <- 
  db %>% 
  mutate(Year = Year %>% as.double()) %>% 
  filter(Country %in% cts20,
         Cause %in% c("1000", "AAA"),
         Year >= 2015) %>% 
  select(-Admin1, -SubDiv, -Cause, -List) %>% 
  arrange(name, Year)

unique(db20$Frmat)

db20_2 <- 
  db20 %>% 
  gather(contains("Deaths"), key = Age, value = Deaths) %>% 
  filter(Age %in% paste0("Deaths", 2:25)) %>% 
  mutate(Age = recode(Age,
                      'Deaths2' = 0,
                      'Deaths3' = 1,
                      'Deaths4' = 1,
                      'Deaths5' = 1,
                      'Deaths6' = 1,
                      'Deaths7' = 5,
                      'Deaths8' = 10,
                      'Deaths9' = 15,
                      'Deaths10' = 20,
                      'Deaths11' = 25,
                      'Deaths12' = 30,
                      'Deaths13' = 35,
                      'Deaths14' = 40,
                      'Deaths15' = 45,
                      'Deaths16' = 50,
                      'Deaths17' = 55,
                      'Deaths18' = 60,
                      'Deaths19' = 65,
                      'Deaths20' = 70,
                      'Deaths21' = 75,
                      'Deaths22' = 80,
                      'Deaths23' = 85,
                      'Deaths24' = 90,
                      'Deaths25' = 95),
         Deaths = Deaths %>% as.double(),
         Country = Country %>% as.double(),
         Sex = recode(Sex,
                      "1" = "m",
                      "2" = "f")) %>%
  drop_na(Deaths) %>% 
  rename(Code = Country,
         Country = name) %>% 
  mutate(Code = countrycode(Country, origin = 'country.name', destination = 'iso3c'),
         Country = recode(Country,
                          "United Kingdom, England and Wales" = "England and Wales",
                          "United Kingdom, Scotland" = "Scotland",
                          "Czech Republic" = "Czechia"),
         Code = case_when(Country == "England and Wales" ~ "GBR-ENW",
                          Country == "Scotland" ~ "GBR-SCO",
                          Country == "Northern Ireland" ~ "GBR-NIR",
                          TRUE ~ Code)) %>% 
  select(-Frmat, -IM_Frmat) %>% 
  group_by(Country, Code, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(age_up = ifelse(Age == max(Age), 110 - Age, lead(Age) - 1),
         Source = "who_mort_db")

unique(db20_2$Code)

write_csv(db20_2, "Output/who.csv")
