library(readxl)
library(tidyverse)
library(countrycode)
library(lubridate)

write_excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,file = paste0("clipboard-", object.size(x)),sep="\t",row.names=row.names,col.names=col.names,...)
}

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

# chunk <- 
#   all_in2 %>% 
#   filter(Source == "unpd",
#          Country == "Belgium", 
#          Year == 2018, 
#          Sex == "f")

# chunk <- temp5

harmon_age_intervals <- function(chunk){
  
  int <- 
    ref_ages %>% 
    semi_join(chunk, by = c("Source", "Country")) %>% 
    pull(Age)
  
  if(max(int) <= 110){
    int <- c(int, 120)
  }
  
  labs <- int[1:length(int)-1]
  chunk %>% 
    mutate(Age = cut(Age, breaks = int, include.lowest = TRUE, right = FALSE, labels = labs),
           Age = as.numeric(as.character(Age))) %>% 
    group_by(Age) %>% 
    summarise(Deaths = sum(Deaths)) %>% 
    ungroup()
}

assign_age_intervals <- function(chunk){
  ct <- unique(chunk$Country)
  int <- 
    ref_ages %>% 
    filter(Country == ct) %>% 
    pull(Age)
  
  if(max(int) <= 110){
    int <- c(int, 110)
  }
  
  labs <- int[1:length(int)-1]
  chunk %>% 
    mutate(Age = cut(Age, breaks = int, include.lowest = TRUE, right = FALSE, labels = labs),
           Age = as.numeric(as.character(Age))) %>% 
    group_by(Age) %>% 
    summarise(Deaths = sum(Deaths)) %>% 
    ungroup()
}


# Groupping population using the same age intervals as mortality data
assign_age_invals_pop <- function(chunk){
  ct <- unique(chunk$Country)
  yr <- unique(chunk$Year)
  sx <- unique(chunk$Sex)
  
  int <- 
    ref_ages %>% 
    filter(Country == ct,
           Year == yr) %>% 
    pull(Age) %>% 
    sort
  
  if(max(int) <= 110){
    int <- c(int, 110)
  }
  
  labs <- int[1:length(int)-1]
  chunk %>% 
    mutate(Age = cut(Age, breaks = int, include.lowest = TRUE, right = FALSE, labels = labs),
           Age = as.numeric(as.character(Age))) %>% 
    group_by(Age) %>% 
    summarise(Population = sum(Population)) %>% 
    ungroup()
}

sum_source <- function(db){
 # test <-
  db %>% 
    group_by(Source, Country, Year, Sex) %>% 
    mutate(ages = n(),
           infd = ifelse(any(Age == 0 & age_spn == 1), 1, 0)) %>% 
    ungroup() %>% 
    group_by(Source, Country, Year, Age) %>% 
    mutate(sexs = n()) %>% 
    ungroup() %>% 
    group_by(Source, Country, Sex, Age) %>% 
    mutate(years = n()) %>% 
    ungroup() %>% 
    group_by(Source, Country, Sex) %>% 
    mutate(years = ifelse(Year >= 2020, max(years), years)) %>% 
    ungroup() %>% 
    group_by(Source, Country) %>% 
    filter(!(sexs == 3 & Sex == "t")) %>% 
    summarise(Deaths = sum(Deaths),
              infd = min(infd),
              ages = min(ages),
              sexs = min(sexs),
              years = min(years),
              period = paste(min(Year), max(Year), sep = "-")) %>% 
    ungroup() %>% 
    unique() %>% 
    group_by(Country) %>% 
    mutate(Sources = n()) %>% 
    select(Country, Sources, everything())
}


