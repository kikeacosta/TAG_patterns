rm(list=ls())
source("R/00_functions.R")

# Loading STMF data
# ~~~~~~~~~~~~~~~~~

# downloading the last version of STMF Mortality input data zip 
# this version as of 25 May 2021
download.file("https://www.mortality.org/File/GetDocument/Public/STMF/Inputs/STMFinput.zip", 
              "Data/STMFinput.zip")

# list of country codes in STMF
zipdf <- unzip("Data/STMFinput.zip", list = TRUE)

# loading all cause deaths from all countries in STMF
db_d <- tibble()
for(i in 1:length(zipdf$Name)){
  csv_file <- zipdf$Name[i]
  print(csv_file)
  temp <- 
    read_csv(unz("Data/STMFinput.zip", csv_file)) %>% 
    mutate(Week = Week %>% as.double,
           Deaths = Deaths %>% as.double())
  db_d <- 
    db_d %>% 
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

# countries with full 2022
cts_2022 <- 
  db_d %>% 
  filter(Year == 2022) %>% 
  group_by(PopCode) %>% 
  filter(max(Week) == 52) %>% 
  pull(PopCode) %>% unique()

# filtering periods 2015-2021 with full annual information
dts <- 
  db_d %>% 
  mutate(PopCode = ifelse(PopCode == "a", "NOR", PopCode)) %>% 
  filter(Year %in% 2010:2021) %>% 
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

# manual homogenization ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# France: leave it as it is
# USA: exclude it
# Scotland: exclude it
# England and Wales: exclude it
# Northern Ireland: exclude it
# Australia: Exclude it because it is hospital data and age groups are horrible
exc <- c("USA", "GBR_NIR", "GBRTENW", "GBR_SCO", "AUS")
unique(dts$PopCode)
dts2 <- 
  dts %>% 
  filter(!PopCode %in% exc) %>% 
  drop_na(Deaths) %>% 
  mutate(Age = case_when(PopCode %in% c("FRATNP") & Age == "1" ~ "0", 
                         PopCode %in% c("FRATNP") & Age == "95" ~ "90", 
                         TRUE ~ Age)) %>% 
  group_by(PopCode, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  filter(Age != "UNK") %>% 
  mutate(Sex = recode(Sex,
                      "b" = "t"),
         PopCode = recode(PopCode,
                          "DEUTNP" = "DEU",
                          "FRATNP" = "FRA",
                          "NZL_NP" = "NZL")) %>%
  rename(Code = PopCode) 

# re-scaling ages and sexes
dts3 <- 
  dts2 %>% 
  group_by(Code, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Code, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Code, Year, Sex, Age)

dts4 <- 
  dts3 %>% 
  mutate(Source = "stmf",
         Country = countrycode(sourcevar = Code, 
                               origin = "iso3c", 
                               destination = "country.name")) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

# saving data
write_csv(dts4, "data_inter/stmf.csv")

dts4 <- read_csv("data_inter/stmf.csv") 

