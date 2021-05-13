
source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")

library(googlesheets4)
library(tidyverse)
library(cartography)
library(rgdal)
library(tmap)
library(sf)
data(World)

# checking coordinate system
st_crs(World)

World$name <- as.character(World$name)

World <- 
  World %>% 
  mutate(name = case_when(name == "Czech Rep." ~ "Czechia",
                          name == "Swaziland" ~ "Eswatini",
                          name == "United States" ~ "USA",
                          TRUE ~ name))
# remove Antarctica
World <- World[World$name != "Antarctica",]

