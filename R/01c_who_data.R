library(readxl)
source("R/00_functions.R")

db <- 
  read_csv("Data/icd_raw.csv",
           col_types = cols(.default = "c")) %>% 
  select(-Country) %>% 
  rename(Country = name)

# selecting: 
# countries with data in 2020; data by age; all-cause mortality; data since 2015
db20 <- 
  db %>% 
  mutate(Year = Year %>% as.double()) %>% 
  group_by(Country) %>% 
  filter(max(Year) >= 2020 &  Frmat != "09") %>% 
  ungroup() %>% 
  filter(Cause %in% c("1000", "AAA"),
         Year >= 2015) %>% 
  select(-Admin1, -SubDiv, -Cause, -List) %>% 
  arrange(Country, Year)

unique(db20$Frmat)

db20_2 <- 
  db20 %>% 
  select(Country, Year, Sex, starts_with("Deaths")) %>% 
  gather(starts_with("Deaths"), key = Age, value = Deaths) %>% 
  filter(Age %in% paste0("Deaths", 2:25)) %>% 
  mutate(Age = recode(Age,
                      'Deaths2' = "TOT",
                      'Deaths2' = "0",
                      'Deaths3' = "1",
                      'Deaths4' = "1",
                      'Deaths5' = "1",
                      'Deaths6' = "1",
                      'Deaths7' = "5",
                      'Deaths8' = "10",
                      'Deaths9' = "15",
                      'Deaths10' = "20",
                      'Deaths11' = "25",
                      'Deaths12' = "30",
                      'Deaths13' = "35",
                      'Deaths14' = "40",
                      'Deaths15' = "45",
                      'Deaths16' = "50",
                      'Deaths17' = "55",
                      'Deaths18' = "60",
                      'Deaths19' = "65",
                      'Deaths20' = "70",
                      'Deaths21' = "75",
                      'Deaths22' = "80",
                      'Deaths23' = "85",
                      'Deaths24' = "90",
                      'Deaths25' = "95"),
         Deaths = Deaths %>% as.double(),
         Sex = recode(Sex,
                      "1" = "m",
                      "2" = "f")) %>%
  drop_na(Deaths) %>% 
  mutate(Code = countrycode(Country, origin = 'country.name', destination = 'iso3c'),
         Country = recode(Country,
                          "United Kingdom, England and Wales" = "England and Wales",
                          "United Kingdom, Scotland" = "Scotland",
                          "Czech Republic" = "Czechia"),
         Code = case_when(Country == "England and Wales" ~ "GBR-ENW",
                          Country == "Scotland" ~ "GBR-SCO",
                          Country == "Northern Ireland" ~ "GBR-NIR",
                          TRUE ~ Code)) %>% 
  group_by(Country, Code, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(Source = "who_mort_db")

unique(db20_2$Code)
unique(db20_2$Age)
unique(db20_2$Sex)

# adding total sex
db20_3 <- 
  db20_2 %>% 
  group_by(Country, Code, Year, Age, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  mutate(Sex = "t") %>% 
  bind_rows(db20_2) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Country, Code, Year, Sex, suppressWarnings(as.integer(Age)))

write_csv(db20_3, "Output/who.csv")
