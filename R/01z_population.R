source("R/00_functions.R")

# Data from the WPP 2022
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/2_Population/WPP2022_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/2_Population/WPP2022_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx

pop_f1 <- 
  read_xlsx("Data/WPP2022_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx",
            skip = 16) %>% 
  select(3, 11:112) %>% 
  rename(Country = 1, 
         Year = 2) %>% 
  filter(Year >= 2010) %>% 
  gather(-Country, -Year, key = "Age", value = "Population") %>% 
  mutate(Sex = "f") %>% 
  arrange(Country, Age)

pop_f2 <- 
  read_xlsx("Data/WPP2022_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE.xlsx",
            sheet = 2,
            skip = 16) %>% 
  select(3, 11:112) %>% 
  rename(Country = 1, 
         Year = 2) %>% 
  filter(Year <= 2025) %>% 
  gather(-Country, -Year, key = "Age", value = "Population") %>% 
  mutate(Sex = "f") %>% 
  arrange(Country, Age)

pop_m1 <- 
  read_xlsx("Data/WPP2022_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx",
            skip = 16) %>% 
  select(3, 11:112) %>% 
  rename(Country = 1, 
         Year = 2) %>% 
  filter(Year >= 2010) %>% 
  gather(-Country, -Year, key = "Age", value = "Population") %>% 
  mutate(Sex = "m") %>% 
  arrange(Country, Age)

pop_m2 <- 
  read_xlsx("Data/WPP2022_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx",
            sheet = 2,
            skip = 16) %>% 
  select(3, 11:112) %>% 
  rename(Country = 1, 
         Year = 2) %>% 
  filter(Year <= 2025) %>% 
  gather(-Country, -Year, key = "Age", value = "Population") %>% 
  mutate(Sex = "m") %>% 
  arrange(Country, Age)

pop <- 
  bind_rows(pop_f1, pop_f2, pop_m1, pop_m2) %>% 
  mutate(Country = case_when(Country == "United States of America" ~ "USA",
                             Country == "Republic of Korea" ~ "South Korea",
                             Country == "Russian Federation" ~ "Russia",
                             Country == "China, Taiwan Province of China" ~ "Taiwan",
                             Country == "China, Hong Kong SAR" ~ "Hong Kong",
                             Country == "Iran (Islamic Republic of)" ~ "Iran",
                             Country == "Republic of Moldova" ~ "Moldova",
                             Country == "Bolivia (Plurinational State of)" ~ "Bolivia",
                             TRUE ~ Country),
         Age = ifelse(Age == "100+", "100", Age),
         Age = as.integer(Age),
         Population = as.double(Population) * 1000) %>% 
  arrange(Country, Year, Sex, Age)

pop2 <- 
  pop %>% 
  group_by(Country, Year, Age) %>% 
  summarise(Population = sum(Population)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>% 
  bind_rows(pop) %>% 
  select(Country, Year, Sex, Age, Population) %>% 
  arrange(Country, Year, Sex, Age)

write_csv(pop2, file = "data_inter/offsets.csv")

