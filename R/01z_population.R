library(tidyverse)
library(readr)
WPP_hist <- read_csv("https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_PopulationBySingleAgeSex_1950-2019.csv")
WPPproj <- read_csv("https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_PopulationBySingleAgeSex_2020-2100.csv")
  # WPP <- read_csv("Data/WPP2019_PopulationBySingleAgeSex_1950-2019.csv") 
offsets <-
WPP %>% 
  pivot_longer(PopMale:PopTotal, names_to = "Sex", values_to = "Population") %>% 
  mutate(Sex = recode(Sex, 
                      "PopMale" = "m",
                      "PopFemale" = "f",
                      "PopTotal" = "b"),
         Region = "All") %>% 
  select(Age = AgeGrpStart, Population, Country = Location, Region,  Sex) %>% 
  mutate(Population = Population * 1000) %>% 
  bind_rows(offsets) %>% 
  arrange(Country, Region, Sex, Age)