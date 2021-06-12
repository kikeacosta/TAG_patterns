library(tidyverse)
library(readxl)

read_excel("Data/south_africa_input.xlsx") %>% 
  pivot_longer(`2015_F`:`2020_M`, names_to = "ys", values_to = "Deaths") %>% 
  separate(ys, into = c("Year","Sex"),sep="_", convert = TRUE) %>% 
  mutate(Sex = tolower(Sex)) %>% 
  pivot_wider(names_from = "Sex", values_from = "Deaths") %>% 
  mutate(t = f + m) %>% 
  pivot_longer(f:t, names_to = "Sex", values_to = "Deaths") %>% 
  mutate(Country = "South Africa",
         Code = "ZAF",
         Source = "SAMRC/UCT") %>% 
  select(Country,Code,Year,Sex,Age,Deaths,Source) %>% 
  write_csv("Output/south_africa.csv")
