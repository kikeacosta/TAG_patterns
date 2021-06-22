library(here)
library(readxl)
library(tidyverse)
library(countrycode)
library(lubridate)


rescale_age <- function(chunk){
  TOT <- chunk %>% dplyr::filter(Age == "TOT") %>% dplyr::pull(Deaths)
  chunk %>% 
    dplyr::filter(Age != "TOT") %>% 
    mutate(Deaths = (Deaths / sum(Deaths)) * TOT)
}

# rescale deaths by sex to total sexes
rescale_sex <- function(chunk){
  TOT <- chunk %>% dplyr::filter(Sex == "t") %>% dplyr::pull(Deaths)
  temp1 <- 
    chunk %>% 
    dplyr::filter(Sex != "t") %>% 
    mutate(Deaths = (Deaths / sum(Deaths)) * TOT) %>% 
    bind_rows(chunk %>% 
                dplyr::filter(Sex == "t"))
  
}

std_db <- function(db){
  db2 <- db %>% 
    select(Country, Year = YearOccurrence, Sex, Age = AgeStart, Deaths)
}
