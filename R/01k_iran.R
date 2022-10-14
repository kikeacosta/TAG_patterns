library(here)
library(tidyverse)
library(countrycode)

iran <- read_csv("Data/iran_deaths_iso_weeks_sex_age.csv")

unique(iran$country_name)
iran2 <- 
  iran %>% 
  separate(age_group, c("Age", "trash"), sep = "-") %>% 
  mutate(Age = case_when(Age == "85+" ~ "85",
                         TRUE ~ Age),
         Sex = str_sub(sex, 1, 1)) %>% 
  group_by(country_name, year, Age, Sex) %>% 
  summarise(Deaths = sum(deaths),
            Code = "IRN",
            .groups = "drop") %>% 
  select(Country = country_name, Code, Year = year, Sex, Age, Deaths)

iran_sex_t <- 
  iran2 %>% 
  group_by(Country, Code, Year, Age) %>% 
  summarise(Deaths = sum(Deaths),
            .groups = "drop") %>% 
  mutate(Sex = "t")

iran_age_t <- 
  iran2 %>% 
  group_by(Country, Code, Year, Sex) %>% 
  summarise(Deaths = sum(Deaths),
            .groups = "drop") %>% 
  mutate(Age = "TOT")

iran_age_sex_t <- 
  iran_age_t %>% 
  group_by(Country, Code, Year, Age) %>% 
  summarise(Deaths = sum(Deaths),
            .groups = "drop") %>% 
  mutate(Sex = "t")

iran_out <- 
  bind_rows(iran2,
            iran_sex_t,
            iran_age_t,
            iran_age_sex_t) %>% 
  filter(Year <= 2021) %>% 
  arrange(Country, Code, Year, Sex, suppressWarnings(as.numeric(Age))) %>% 
  mutate(Source = "wmd") %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Country, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Code, Year, Sex, Age)

write_csv(iran_out, "Output/iran.csv")
