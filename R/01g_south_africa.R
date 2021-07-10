library(tidyverse)
library(readxl)

zaf <- read_excel("Data/south_africa_input.xlsx") %>% 
  pivot_longer(`2015_F`:`2020_M`, names_to = "ys", values_to = "Deaths") %>% 
  separate(ys, into = c("Year","Sex"),sep="_", convert = TRUE) %>% 
  mutate(Sex = tolower(Sex)) %>% 
  pivot_wider(names_from = "Sex", values_from = "Deaths") %>% 
  mutate(t = f + m) %>% 
  pivot_longer(f:t, names_to = "Sex", values_to = "Deaths") %>% 
  mutate(Country = "South Africa",
         Code = "ZAF",
         Source = "SAMRC/UCT") %>% 
  select(Country,Code,Year,Sex,Age,Deaths,Source)

zaf2 <- 
  zaf %>% 
  tidyr::complete(Country, Code, Year, Sex, Age, fill = list(Deaths = 0)) %>% 
  group_by(Country, Code, Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(zaf %>% 
              mutate(Age = as.character(Age))) %>% 
  arrange(Year, Sex) %>% 
  mutate(Source = "samrc/uct")

write_csv(zaf2, "Output/south_africa.csv")
