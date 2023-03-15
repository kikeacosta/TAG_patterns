rm(list=ls())
source("R/00_functions.R")

# loading pre-processed data
usa1  <- 
  read_tsv("Data/USA/wonder_10_20.txt") %>% 
  select(Age = 3,
         Sex = 5,
         Year = 6,
         Deaths) %>% 
  drop_na(Year) %>% 
  mutate(Age = Age %>% as.character())


usa2 <- 
  read_tsv("Data/USA/wonder_18_21.txt") %>%
  select(Age = 3,
         Sex = 5,
         Year = 6,
         Deaths) %>% 
  filter(Year == 2021) %>% 
  # drop_na(Age) %>% 
  mutate(Age = Age %>% as.character())

unique(usa1$Sex)

usa <- 
  bind_rows(usa1, usa2)

usa_tot <- 
  usa %>% 
  group_by(Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(usa %>% filter(!is.na(Age) & Age !="NS")) %>% 
  group_by(Year, Sex) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double())

usa_out <- 
  usa_tot %>% 
  group_by(Year, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>% 
  bind_rows(usa_tot) %>% 
  mutate(Sex = str_to_lower(Sex),
         Source = "direct",
         Country = "USA",
         Code = "USA") %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

write_csv(usa_out, "data_inter/usa.csv")

