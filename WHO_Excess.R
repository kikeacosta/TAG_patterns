library(readxl)
library(tidyverse)
library(countrycode)
WHOin <- read_excel("Data/WHO_Allcause_Mortality_Data_Call_10.05.2021.xlsx")

str(WHOin)
iso3 <- WHOin$country %>% unique()

WHOin %>% 
  filter(country == "ZMB",
         sex == "Total") %>% 
  group_by(year) %>% 
  summarize(d = sum(deaths))
WHOin %>% 
  filter(country == "ZMB",
         sex == "Total",
         time_unit !="annual") %>% 
  group_by(year) %>% 
  summarize(d = sum(deaths))
